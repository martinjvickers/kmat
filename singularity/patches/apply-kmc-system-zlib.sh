#!/usr/bin/env python3
"""Patch refresh-bio/KMC Makefile (v3.2.4) to link system zlib instead of Cloudflare."""
from __future__ import annotations

from pathlib import Path

mk = Path("Makefile")
text = mk.read_text()

replacements = [
    (
        "STATIC_CFLAGS = -static -Wl,--whole-archive -lpthread -Wl,--no-whole-archive",
        "STATIC_CFLAGS = -static-libgcc -static-libstdc++ -pthread",
    ),
    (
        "STATIC_LFLAGS = -static -Wl,--whole-archive -lpthread -Wl,--no-whole-archive",
        "STATIC_LFLAGS = -static-libgcc -static-libstdc++ -pthread -lz",
    ),
    ("LIB_ZLIB=3rd_party/cloudflare/libz.a", "LIB_ZLIB="),
    (
        "$(CC) $(CFLAGS) -I 3rd_party/cloudflare -c $< -o $@",
        "$(CC) $(CFLAGS) -c $< -o $@",
    ),
    (
        "kmc: $(KMC_CLI_OBJS) $(LIB_KMC_CORE) $(LIB_ZLIB)",
        "kmc: $(KMC_CLI_OBJS) $(LIB_KMC_CORE)",
    ),
    (
        "kmc_tools: $(KMC_TOOLS_OBJS) $(KMC_API_OBJS) $(KFF_OBJS) $(LIB_ZLIB)",
        "kmc_tools: $(KMC_TOOLS_OBJS) $(KMC_API_OBJS) $(KFF_OBJS)",
    ),
    (
        "$(CC) $(CLINK) -I 3rd_party/cloudflare -o $(OUT_BIN_DIR)/$@ $^",
        "$(CC) $(CLINK) -o $(OUT_BIN_DIR)/$@ $^",
    ),
    (
        "dummy := $(shell git submodule update --init --recursive)",
        "dummy :=",
    ),
]

for old, new in replacements:
    if old not in text:
        raise SystemExit(f"Makefile patch failed: expected text not found:\n{old}")
    text = text.replace(old, new, 1)

mk.write_text(text)
print("patched Makefile for system zlib")
