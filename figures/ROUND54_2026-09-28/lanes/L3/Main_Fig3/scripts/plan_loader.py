"""Round 54, lane L3: import the planner's PLAN (00_plan_r54.py, a module name that starts with a digit) without running its report."""
import importlib, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
PLAN = importlib.import_module("00_plan_r54").PLAN
