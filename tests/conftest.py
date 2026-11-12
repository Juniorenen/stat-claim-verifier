"""Windows compatibility patch for genlayer-test Direct Mode.

genlayer-test 0.29.2 injects the transaction message by replacing fd 0 with a
temp file, then calls os.unlink() on that path. POSIX allows unlinking an open
file; Windows raises WinError 32. The temp file therefore leaks AND the test
aborts. We replace the injector with an equivalent that defers deletion to
interpreter exit.

This patch is test-harness only: it never touches the contract under test.
"""

import atexit
import os
import tempfile

import gltest.direct.loader as _loader


def _inject_message_to_fd0_windows(vm) -> None:
    from genlayer.py import calldata
    from genlayer.py.types import Address

    sender_addr = vm.sender
    if isinstance(sender_addr, bytes):
        sender_addr = Address(sender_addr)

    contract_addr = vm._contract_address
    if isinstance(contract_addr, bytes):
        contract_addr = Address(contract_addr)

    origin_addr = vm.origin
    if isinstance(origin_addr, bytes):
        origin_addr = Address(origin_addr)

    message_data = {
        "contract_address": contract_addr,
        "sender_address": sender_addr,
        "origin_address": origin_addr,
        "stack": [],
        "value": vm._value,
        "datetime": vm._datetime,
        "is_init": False,
        "chain_id": vm._chain_id,
        "entry_kind": 0,
        "entry_data": b"",
        "entry_stage_data": None,
    }

    encoded = calldata.encode(message_data)

    fd, path = tempfile.mkstemp()
    os.write(fd, encoded)
    os.lseek(fd, 0, os.SEEK_SET)

    vm._original_stdin_fd = os.dup(0)
    os.dup2(fd, 0)
    os.close(fd)

    # Windows cannot unlink a file that is still open (fd 0 holds it).
    atexit.register(lambda p=path: os.path.exists(p) and os.remove(p))


if os.name == "nt":
    _loader._inject_message_to_fd0 = _inject_message_to_fd0_windows
