import os
from pathlib import Path

root = Path(".")

# Find BC files
print("BC_P01 files:")
for f in root.rglob("*BC*P01*"):
    if f.is_file() and ("norm" in f.name or "filtered" in f.name):
        print(f"  {f}")

print("\nLUAD_P01 files:")
for f in root.rglob("*LUAD*P01*"):
    if f.is_file() and ("norm" in f.name or "filtered" in f.name):
        print(f"  {f}")

print("\nPRAD_P01 files:")
for f in root.rglob("*PRAD*P01*"):
    if f.is_file() and ("norm" in f.name or "filtered" in f.name):
        print(f"  {f}")
