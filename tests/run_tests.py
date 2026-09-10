"""Runs the theorem-tests of fieldbit (exact arithmetic). Exit code 1 if any fails."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from fieldbit import run_tests
ok_all = True
for name, val, ok in run_tests():
    print(("OK  " if ok else "FAIL"), name, "=", val); ok_all &= bool(ok)
sys.exit(0 if ok_all else 1)
