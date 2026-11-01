"""Assemble a GIF from the real Qt frames written by nucleus_viewer --capture."""
from pathlib import Path
from PIL import Image
root = Path(__file__).resolve().parent.parent
folder = root / "results/screenshots"
frames = [Image.open(folder / "frames" / f"frame_{i}.png").convert("RGB") for i in range(36)]
if any(frame.size != frames[0].size for frame in frames):
    raise ValueError("Frame sizes differ")
frames[0].save(folder / "geant4_rotation.gif", save_all=True, append_images=frames[1:],
               duration=130, loop=0, optimize=True)
print("Saved GIF from 36 captured Geant4 Qt frames")
