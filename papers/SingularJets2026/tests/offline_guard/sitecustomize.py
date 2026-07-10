"""Fail immediately if an offline reproduction process opens a network socket."""

from __future__ import annotations

import socket


def _blocked(*_args, **_kwargs):
    raise RuntimeError("network access is disabled for the offline capsule check")


socket.create_connection = _blocked
socket.socket.connect = _blocked
socket.socket.connect_ex = _blocked
