# Unit tests for Scripts/tools.py

# These run real shell commands (true, exit, sh loops) through tools.run_tool, so no
# bioinformatics programs need to be installed.

import pytest

import tools


@pytest.fixture(autouse=True)
def clean_tools_state():
	tools.reset()
	yield
	tools.reset()


# --------------------------------------------------------------------------------------
# test run_tool: success and failure detection
# --------------------------------------------------------------------------------------

def test_success_returns_true():
	assert tools.run_tool('true', 'X', 'test') is True
	assert tools.FAILURES == []


def test_nonzero_exit_is_recorded(capsys):
	ok = tools.run_tool("sh -c 'echo bad >&2; exit 2'", 'MAFFT', 'guidance', 'OG1')

	assert ok is False
	assert len(tools.FAILURES) == 1
	failure = tools.FAILURES[0]
	assert failure.tool == 'MAFFT'
	assert failure.step == 'guidance'
	assert failure.gene_family == 'OG1'
	assert failure.exit_code == 2
	assert 'bad' in failure.stderr_tail
	assert 'WARNING: MAFFT failed for gene family OG1 (exit code 2)' in capsys.readouterr().out


def test_empty_output_is_a_failure(tmp_path):
	# A shell redirect creates the output file before the program runs, so a failed
	# program can leave an empty file behind
	out = tmp_path / 'aligned.fasta'

	ok = tools.run_tool('true > ' + str(out), 'MAFFT', 'guidance', 'OG1', outputs=[str(out)])

	assert ok is False
	assert tools.FAILURES[0].exit_code == 0
	assert tools.FAILURES[0].missing_outputs == [str(out)]


def test_missing_output_is_a_failure(tmp_path):
	out = tmp_path / 'never_written.fasta'

	ok = tools.run_tool('true', 'MAFFT', 'guidance', 'OG1', outputs=[str(out)])

	assert ok is False
	assert tools.FAILURES[0].missing_outputs == [str(out)]


def test_existing_output_is_a_success(tmp_path):
	out = tmp_path / 'aligned.fasta'

	ok = tools.run_tool('echo ">a" > ' + str(out), 'MAFFT', 'guidance', 'OG1', outputs=[str(out)])

	assert ok is True
	assert tools.FAILURES == []


# --------------------------------------------------------------------------------------
# test run_tool: the program's error output
# --------------------------------------------------------------------------------------

def test_stderr_tail_keeps_last_20_lines():
	cmd = 'i=0; while [ $i -lt 100 ]; do echo $i >&2; i=$((i+1)); done; exit 1'

	tools.run_tool(cmd, 'IQ-TREE', 'trees', 'OG1')

	assert tools.FAILURES[0].stderr_tail.splitlines() == [str(i) for i in range(80, 100)]


def test_stderr_is_echoed(capsys):
	# The program's error output is captured, then passed on so the SLURM log still has it
	tools.run_tool("sh -c 'echo progress >&2'", 'MAFFT', 'guidance', 'OG1')

	assert 'progress' in capsys.readouterr().err


def test_invalid_utf8_stderr():
	ok = tools.run_tool("printf '\\377\\376 x' >&2; exit 1", 'MAFFT', 'guidance', 'OG1')

	assert ok is False
	assert '�' in tools.FAILURES[0].stderr_tail


# --------------------------------------------------------------------------------------
# test run_tool: when the whole run stops
# --------------------------------------------------------------------------------------

def test_command_not_found_stops(capsys):
	# Exit code 127: the program isn't installed or its module isn't loaded, so every
	# gene family would fail the same way
	with pytest.raises(SystemExit) as excinfo:
		tools.run_tool('no_such_tool_xyz', 'MAFFT', 'guidance', 'OG1')

	assert excinfo.value.code == 1
	out = capsys.readouterr().out
	assert 'ERROR: MAFFT failed for gene family OG1 (exit code 127)' in out
	assert len(tools.FAILURES) == 1


def test_fail_fast_stops():
	tools.configure(True)

	with pytest.raises(SystemExit) as excinfo:
		tools.run_tool('exit 1', 'MAFFT', 'guidance', 'OG1')

	assert excinfo.value.code == 1


def test_critical_stops():
	with pytest.raises(SystemExit) as excinfo:
		tools.run_tool('exit 1', 'cp', 'guidance', critical=True)

	assert excinfo.value.code == 1


def test_failure_without_gene_family_names_the_step(capsys):
	tools.run_tool('exit 1', 'cp', 'guidance')

	assert 'WARNING: cp failed during guidance (exit code 1)' in capsys.readouterr().out


# --------------------------------------------------------------------------------------
# test the failure report
# --------------------------------------------------------------------------------------

def test_write_failure_report(tmp_path):
	tools.run_tool("echo first >&2; echo second >&2; exit 1", 'MAFFT', 'guidance', 'OG1')
	tools.run_tool('exit 3', 'IQ-TREE', 'trees', 'OG2', outputs=[str(tmp_path / 'OG2.treefile')])

	path = tools.write_failure_report(str(tmp_path))

	assert path == str(tmp_path / 'Output' / 'FailedSteps.tsv')
	lines = open(path).read().splitlines()
	assert lines[0] == 'step\tgene_family\ttool\texit_code\tmissing_outputs\tmessage'
	assert len(lines) == 3
	# One line per failure: the message column holds the last line of the error output
	assert lines[1].split('\t') == ['guidance', 'OG1', 'MAFFT', '1', '', 'second']
	assert lines[2].split('\t') == ['trees', 'OG2', 'IQ-TREE', '3', str(tmp_path / 'OG2.treefile'), '']


def test_write_failure_report_without_failures(tmp_path):
	assert tools.write_failure_report(str(tmp_path)) is None
	assert not (tmp_path / 'Output').exists()


def test_failed_gene_families():
	tools.run_tool('exit 1', 'MAFFT', 'guidance', 'OG1')
	tools.run_tool('exit 1', 'trimAl', 'guidance', 'OG1')
	tools.run_tool('exit 1', 'IQ-TREE', 'trees', 'OG2')
	tools.run_tool('exit 1', 'cp', 'guidance')

	assert tools.failed_gene_families() == {'OG1', 'OG2'}


# --------------------------------------------------------------------------------------
# test the fake_shell fixture used by the per-script tests
# --------------------------------------------------------------------------------------

def test_fake_shell(fake_shell, tmp_path):
	(tmp_path / 'a.txt').write_text('x')
	fake_shell.fail('mafft', exit_code=2, stderr='out of memory')

	assert tools.run_tool('mafft in.fasta > out.fasta', 'MAFFT', 'guidance', 'OG1') is False
	assert tools.FAILURES[0].exit_code == 2
	assert 'out of memory' in tools.FAILURES[0].stderr_tail
	# Other programs "succeed" without running; cp/mv/rm run for real
	assert tools.run_tool('iqtree2 -s x', 'IQ-TREE', 'trees', 'OG1') is True
	assert tools.run_tool('cp ' + str(tmp_path / 'a.txt') + ' ' + str(tmp_path / 'b.txt'), 'cp', 'trees') is True
	assert (tmp_path / 'b.txt').read_text() == 'x'
	assert len(fake_shell.calls) == 3
