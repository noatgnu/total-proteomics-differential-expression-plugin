#!/usr/bin/env python3
"""Merges each platform's deps-partial.json (downloaded into deps-downloaded/deps-<platform>/) into one dependencies.json."""

import glob
import json

merged = {"generatedAt": None, "packages": [], "systemRequirements": {}}

for path in sorted(glob.glob("deps-downloaded/deps-*/deps-partial.json")):
    platform = path.split("deps-downloaded/deps-")[1].split("/")[0]
    with open(path) as f:
        data = json.load(f)
    if data.get("generatedAt"):
        merged["generatedAt"] = data["generatedAt"]
    if not merged["packages"] and data.get("packages"):
        merged["packages"] = data["packages"]
    merged["systemRequirements"][platform] = data.get("systemRequirements", [])

with open("dependencies.json", "w") as f:
    json.dump(merged, f, indent=2)
    f.write("\n")

print(f"Wrote dependencies.json with platforms: {list(merged['systemRequirements'].keys())}")
