from pathlib import Path
import subprocess
import sys
import contextlib


@contextlib.contextmanager
def noop():
    """
    Basically we just need a None 'context' for with() statements
    Thanks StackOverflow
    https://stackoverflow.com/questions/28680442/python-context-for-file-or-none
    """
    yield None


def get_tool_path(path_string, tool_name):
    tool_path = Path(path_string)
    if tool_path.name != tool_name:
        tool_path += tool_name

    if tool_path.exists():
        return tool_path
    else:
        raise FileNotFoundError


def subprocess_wrapper(command_list, stdout_file, stderr_file=None):
    """
    Error and output handling for subprocess.run()
    Takes a list of strings for command and options, and strings/paths for stdout, and optionally stderr.
    If stderr_file or stdout_file is None, it will not be captured
    """
    stderr_location = subprocess.PIPE if stderr_file is not None else None
    try:
        with open(stdout_file, "w") if stdout_file is not None else noop() as output:
            process_result = subprocess.run(command_list, stdout=output, stderr=stderr_location, check=True)
        if stderr_file is not None:
            with open(stderr_file, "w") as error_file:
                error_file.write(str(process_result.stderr))
    except subprocess.CalledProcessError as e:
        print(f"Could not run command :\n {command_list}", file=sys.stderr)
        print(e, file=sys.stderr)
        sys.exit(1)
    except FileNotFoundError as e:
        print(f"Could not find/create output files for command:\n {command_list}", file=sys.stderr)
        print(e, file=sys.stderr)
        sys.exit(1)
