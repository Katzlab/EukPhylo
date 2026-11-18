import os, subprocess, sys
from unittest import mock

import pytest

# The pipeline scripts import each other as top-level modules (e.g. `import guidance`),
# so Scripts/ has to be on the path for the tests to import them.
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'Scripts'))

import tools


class FakeShell:
	"""Stand-in for subprocess.run, so tests can make an external program fail without installing it.

	Commands containing a substring registered with fail() return that exit code and error
	output. Other commands whose first word is cp, mv or rm run for real (so input files
	still get copied); everything else "succeeds" without running. Every command is
	recorded in `calls`.
	"""

	def __init__(self, real_run):
		self.real_run = real_run
		self.failures = {}
		self.calls = []

	def fail(self, substring, exit_code = 1, stderr = 'boom'):
		"""Make every command containing `substring` fail with this exit code and error output."""
		self.failures[substring] = (exit_code, stderr)

	def __call__(self, cmd, *args, **kwargs):
		self.calls.append(cmd)
		for substring, (exit_code, stderr) in self.failures.items():
			if substring in cmd:
				return subprocess.CompletedProcess(cmd, exit_code, stderr = stderr.encode())
		if cmd.split()[0] in ('cp', 'mv', 'rm'):
			return self.real_run(cmd, *args, **kwargs)
		return subprocess.CompletedProcess(cmd, 0, stderr = b'')


@pytest.fixture
def fake_shell():
	"""Route tools.run_tool through a FakeShell for one test, with no failures recorded before or after."""
	tools.reset()
	shell = FakeShell(subprocess.run)
	with mock.patch('tools.subprocess.run', side_effect = shell):
		yield shell
	tools.reset()
