"""Restricted data grammar. No downloaded program is ever executed.

GAP: 'return' nested comma lists of signed decimal integers, then ';'.
PARI: an optional single identifier assignment, then lists/matrices, then ';'.
Semicolons inside brackets separate matrix rows. A backslash-newline is a
PARI physical line continuation. No expressions, calls, comments or holes.
"""

import re

from .exact import Invalid, require


class Parser:
    def __init__(self, text):
        require(len(text) <= 8_000_000, "source exceeds 8 MB limit")
        self.text = text
        self.pos = 0

    def ws(self):
        while self.pos < len(self.text) and self.text[self.pos].isspace():
            self.pos += 1

    def take(self, token):
        self.ws()
        require(self.text.startswith(token, self.pos), f"expected {token!r} at {self.pos}")
        self.pos += len(token)

    def value(self, depth=0):
        require(depth <= 5, "excessive nesting")
        self.ws()
        if self.text.startswith("[", self.pos):
            self.pos += 1
            self.ws()
            if self.text.startswith("]", self.pos):
                self.pos += 1
                return []
            row = [self.value(depth + 1)]
            rows = []
            while True:
                self.ws()
                require(self.pos < len(self.text), "unterminated list")
                sep = self.text[self.pos]
                self.pos += 1
                if sep == "]":
                    return rows + [row] if rows else row
                require(sep in ",;", f"unsupported syntax at {self.pos - 1}")
                if sep == ";":
                    require(all(type(x) is int for x in row), "matrix row must contain integers")
                    rows.append(row)
                    row = []
                row.append(self.value(depth + 1))
        m = re.match(r"[+-]?[0-9]+", self.text[self.pos:])
        require(m is not None and len(m[0]) <= 100, f"expected integer at {self.pos}")
        self.pos += len(m[0])
        return int(m[0])


def parse_source(text, fmt):
    require(fmt in ("gap", "pari", "pari-records"), "unknown source format")
    if fmt == "pari-records":
        return parse_pari_records(text)
    if fmt == "pari":
        text = text.replace("\\\r\n", "").replace("\\\n", "")
    p = Parser(text)
    if fmt == "gap":
        p.take("return")
    else:
        p.ws()
        m = re.match(r"[A-Za-z_][A-Za-z_0-9]*\s*=", text[p.pos:])
        if m:
            p.pos += len(m[0])
    result = p.value()
    p.take(";")
    p.ws()
    require(p.pos == len(text), f"trailing source input at {p.pos}")
    # GAP admits only nested lists, never matrix semicolons.
    if fmt == "gap":
        require(";" not in text[:text.rfind(";")], "matrix syntax is not GAP data")
    return result


def parse_pari_records(text):
    """Actual DATA_pari format: consecutive U1=[...], U2=[...], ... .

    Only monotonically numbered assignments of integer lists are accepted.
    A semicolon between assignments is optional. No assignment is executed.
    """
    p = Parser(text)
    result = []
    while True:
        p.ws()
        if p.pos == len(text):
            break
        p.take(f"U{len(result) + 1}")
        p.take("=")
        result.append(p.value())
        p.ws()
        if text.startswith(";", p.pos):
            p.pos += 1
    require(result, "empty PARI assignment list")
    return result
