"""Capture an actual X11 application window without unrelated desktop content."""
import argparse
from pathlib import Path
import struct
import subprocess
import tempfile
import zlib

def capture(title, output):
    with tempfile.TemporaryDirectory() as temp:
        raw = Path(temp) / "window.xwd"
        subprocess.run(["xwd", "-silent", "-name", title, "-out", str(raw)], check=True, timeout=20)
        data = raw.read_bytes()
    h = struct.unpack(">25I", data[:100])
    header, version, fmt, depth, width, height, _, byte_order, _, _, _, bpp, stride, _, red, green, blue, _, _, colors, *_ = h
    if version != 7 or fmt != 2 or bpp != 32 or (red, green, blue) != (0xff0000, 0xff00, 0xff):
        raise ValueError("Unsupported XWD format; refuse to produce a misleading capture")
    start = header + colors*12
    rows = []
    for y in range(height):
        row = bytearray([0])
        for x in range(width):
            offset = start + y*stride + x*4
            pixel = int.from_bytes(data[offset:offset+4], "little" if byte_order == 0 else "big")
            row.extend(((pixel >> 16) & 255, (pixel >> 8) & 255, pixel & 255))
        rows.append(row)
    def chunk(tag, payload):
        return struct.pack(">I", len(payload)) + tag + payload + struct.pack(">I", zlib.crc32(tag+payload) & 0xffffffff)
    png = b"\x89PNG\r\n\x1a\n" + chunk(b"IHDR", struct.pack(">2I5B", width, height, 8, 2, 0, 0, 0))
    png += chunk(b"IDAT", zlib.compress(b"".join(rows))) + chunk(b"IEND", b"")
    Path(output).write_bytes(png)
    print(f"Captured {title}: {width}x{height} -> {output}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("title")
    parser.add_argument("output")
    args = parser.parse_args()
    capture(args.title, args.output)
