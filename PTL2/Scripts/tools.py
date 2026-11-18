"""Run the external programs EukPhylo part 2 depends on and check that each one worked.

A step has failed when the program exits with a non-zero code, or when an output file it should have written is missing or empty. 
Failures are recorded so that eukphylo.py can list them in Output/FailedSteps.tsv and end the run with a non-zero exit code.
"""

#Dependencies
import os, subprocess, sys
from dataclasses import dataclass, field

#Number of lines of a failed program's error output kept in messages and in the report
STDERR_TAIL_LINES = 20

#Exit code the shell uses when a command is not found (program not installed or module not loaded)
COMMAND_NOT_FOUND = 127

REPORT_COLUMNS = ['step', 'gene_family', 'tool', 'exit_code', 'missing_outputs', 'message']


@dataclass
class ToolFailure:
	"""One failed external program call.

	Args:
	    step (str): Pipeline step that ran the program, e.g. 'guidance' or 'trees'.
	    gene_family (str): Gene family being processed, e.g. 'OG6_100457', or '' for whole-run commands.
	    tool (str): Name of the program, e.g. 'MAFFT'.
	    exit_code (int): The program's exit code (0 if it exited normally but an output is missing).
	    missing_outputs (list of str): Expected output files that were missing or empty.
	    stderr_tail (str): The last lines of the program's error output.
	"""

	step: str
	gene_family: str
	tool: str
	exit_code: int
	missing_outputs: list = field(default_factory = list)
	stderr_tail: str = ''

FAILURES = []
FAIL_FAST = False


def configure(fail_fast):
    """Set whether the first external program failure stops the whole run.

    Args:
        fail_fast (bool): Value of the --fail_fast command-line option.
    """
	
    global FAIL_FAST
    FAIL_FAST = fail_fast


def reset():
	"""Forget all recorded failures and turn --fail_fast off (used between tests)."""
	
	global FAIL_FAST
	FAILURES.clear()
	FAIL_FAST = False


def run_tool(cmd, tool, step, gene_family = '', outputs = (), critical = False):
	"""Run an external program through the shell and check that it worked.

	The command runs with /bin/sh, exactly as os.system would. Its error output is
	captured, then written to stderr once the program finishes, so the SLURM log keeps it.

	Args:
	    cmd (str): Shell command, unchanged from the original os.system call.
	    tool (str): Name of the program, used in messages and the report, e.g. 'MAFFT'.
	    step (str): Pipeline step running the program, e.g. 'guidance'.
	    gene_family (str, optional): Gene family being processed; '' for commands that apply to the whole run.
	    outputs (sequence of str, optional): Files the program must write. A missing or empty file counts as a failure.
	    critical (bool, optional): If True, a failure stops the whole run (e.g. copying the input files).

	Returns (bool): True if the program succeeded. False if it failed and the run continues; the caller should then skip the rest of this gene family.

	Raises:
	    (SystemExit) With code 1 if the step failed and `critical` is True, --fail_fast is on, or the
		command was not found (exit code 127).

	Examples:
	>>> run_tool('mafft in.fasta > out.fasta', 'MAFFT', 'guidance', 'OG6_100457', outputs = ['out.fasta'])  # doctest: +SKIP
	True
	>>> run_tool('exit 2', 'MAFFT', 'guidance', 'OG6_100457')  # doctest: +ELLIPSIS
	<BLANKLINE>
	WARNING: MAFFT failed for gene family OG6_100457 (exit code 2). ...
	False
	"""
	
	result = subprocess.run(cmd, shell = True, stderr = subprocess.PIPE)
	stderr = result.stderr.decode(errors = 'replace') if result.stderr else ''
	if stderr != '':
		sys.stderr.write(stderr)

	missing = [path for path in outputs if not os.path.isfile(path) or os.path.getsize(path) == 0]
	if result.returncode == 0 and len(missing) == 0:
		return True

	failure = ToolFailure(step, gene_family, tool, result.returncode, missing, '\n'.join(stderr.splitlines()[-STDERR_TAIL_LINES:]))
	FAILURES.append(failure)

	stop = critical or FAIL_FAST or result.returncode == COMMAND_NOT_FOUND
	print(failure_message(failure, stop))
	if stop:
		sys.exit(1)

	return False


def failure_message(failure, stop):
	"""Build the WARNING (run continues) or ERROR (run stops) message for a failure.

	Arg:
	    failure (ToolFailure): The failure to describe.
	    stop (bool): Whether the run stops because of this failure.

	Returns:
		(str): The message, starting with a blank line like the pipeline's other messages.

	Examples
	>>> message = failure_message(ToolFailure('guidance', 'OG1', 'MAFFT', 2, [], 'bad'), False)
	>>> message.splitlines()[1]
	'WARNING: MAFFT failed for gene family OG1 (exit code 2). Skipping the rest of this gene family; it will be listed in Output/FailedSteps.tsv.'
	>>> message.splitlines()[-1].strip()
	'bad'
	"""
	
	if failure.gene_family != '':
		where = 'for gene family ' + failure.gene_family
		skipped = 'this gene family'
	else:
		where = 'during ' + failure.step
		skipped = 'this step'

	summary = failure.tool + ' failed ' + where + ' (exit code ' + str(failure.exit_code) + ').'
	if stop:
		lines = ['\nERROR: ' + summary + ' Stopping the run.']
		if failure.exit_code == COMMAND_NOT_FOUND:
			lines.append('\tThe command was not found. Check that ' + failure.tool + ' is installed and its module is loaded.')
	else:
		lines = ['\nWARNING: ' + summary + ' Skipping the rest of ' + skipped + '; it will be listed in Output/FailedSteps.tsv.']

	for path in failure.missing_outputs:
		lines.append('\tMissing or empty output: ' + path)
	if failure.stderr_tail != '':
		lines.append('\tLast lines of ' + failure.tool + "'s error output:")
		lines.extend(['\t\t' + line for line in failure.stderr_tail.splitlines()])

	return '\n'.join(lines)


def failed_gene_families():
	"""Return the names of the gene families that had at least one failed step.

	Returns:
	    (set of str): Gene family names; whole-run failures (no gene family) are not included.
	"""
	return { failure.gene_family for failure in FAILURES if failure.gene_family != '' }


def write_failure_report(output_dir):
	"""Write the recorded failures to <output_dir>/Output/FailedSteps.tsv.

	Args:
	    output_dir (str): The --output folder of the run.

	Returns:
	    (str or None): Path of the report, or None if nothing failed (no file is written).

	Examples:
	>>> step	gene_family	tool	exit_code	missing_outputs	message
		guidance	OG6_100457	MAFFT	1		Out of memory
	"""
	if len(FAILURES) == 0:
		return None

	os.makedirs(output_dir + '/Output', exist_ok = True)
	path = output_dir + '/Output/FailedSteps.tsv'
	with open(path, 'w') as o:
		o.write('\t'.join(REPORT_COLUMNS) + '\n')
		for failure in FAILURES:
			stderr_lines = [line.strip() for line in failure.stderr_tail.splitlines() if line.strip() != '']
			message = stderr_lines[-1].replace('\t', ' ') if len(stderr_lines) > 0 else ''
			o.write('\t'.join([failure.step, failure.gene_family, failure.tool, str(failure.exit_code), ','.join(failure.missing_outputs), message]) + '\n')

	return path
