#!/usr/bin/env python3
"""ONE rounding rule for the T90 paper: half up, through Decimal, on the source value.

Canonical copy: ROUND30_2026-09-08/V7_L11_ROUNDING/rounding_v7.py
Importable copy: T90_Manuscript/numbers/rounding_v7.py (numbers/ is already on sys.path in
every figure builder and in the tables lane's lib_tables_v7, so both code paths can
`import rounding_v7`). The two copies are byte for byte the same file.

WHY THIS EXISTS
The figure builders print with Decimal(repr(float(v))).quantize(..., ROUND_HALF_UP).
The table generator prints with an f-string, f"{v:.2f}", which is Python's round-half-to-even
on the binary double. The two rules part on every value that sits on an exact half at the
printed number of decimals: 3.695 prints 3.70 on a figure and 3.69 in a table. Half up is the
rule the paper uses everywhere else (P values, q values, percentages, Supp Table 14) and it is
the rule the figures already follow, so half up is the one rule.

HOW TO USE IT
    from rounding_v7 import half_up, hr_ci, p_ama, p3, pct, strip0
    half_up(3.695, 2)            -> '3.70'
    hr_ci(7.104, 3.695, 13.659)  -> '7.10 (3.70-13.66)'

RULES
1. Apply it to the SOURCE value, never to a value that has already been rounded once.
   Double rounding moves digits: 1.2349 -> 1.235 -> 1.24, while the source rounds to 1.23.
   Every hazard-ratio numbers file of this paper stores round(x, 3), so a two-decimal print
   from those files is a second rounding. `verify/full_precision_check` records where that
   second rounding could change a printed digit.
2. Decimal(repr(float(v))) and not Decimal(v): repr gives the shortest decimal string that
   round-trips the double, which is the number a reader means when they write 3.695. The raw
   binary double for 3.695 is 3.694999999999999951..., so Decimal(v) would round it DOWN and
   silently reintroduce the defect this module removes.
3. Never str(v) on a numpy scalar: numpy repr can carry a dtype tag. float(v) first.
"""
from decimal import Decimal, ROUND_HALF_UP, localcontext

__all__ = ["half_up", "half_up_dec", "hr_ci", "p_ama", "p3", "pct", "strip0", "RULE"]

RULE = "Decimal(repr(float(v))).quantize(1e-dp, ROUND_HALF_UP)"


def half_up_dec(v, dp):
    """Decimal of v rounded half up to dp decimals."""
    if v is None:
        raise ValueError("half_up_dec on None")
    f = float(v)
    if f != f:
        raise ValueError("half_up_dec on NaN")
    if f in (float("inf"), float("-inf")):
        raise ValueError("half_up_dec on an infinite value")
    d = Decimal(repr(f))
    # quantize raises InvalidOperation when the result needs more digits than the context allows,
    # which happens on the huge upper limits of an unstable fit (1e+30). Give it room rather than
    # letting the caller silently fall back to a different rule.
    with localcontext() as ctx:
        ctx.prec = max(60, d.adjusted() + int(dp) + 10)
        return d.quantize(Decimal(1).scaleb(-int(dp)), rounding=ROUND_HALF_UP)


def half_up(v, dp):
    """String of v rounded half up to dp decimals, with dp digits always shown."""
    return str(half_up_dec(v, dp))


def hr_ci(hr, lo, hi, dp=2, sep="-"):
    """A hazard ratio, odds ratio or risk ratio with its interval, all half up."""
    return f"{half_up(hr, dp)} ({half_up(lo, dp)}{sep}{half_up(hi, dp)})"


def strip0(s):
    return s[1:] if s.startswith("0.") else s


def p_ama(v):
    """P or q in the style the supplement already prints (AMA, half up).
    <.001 below .001, three digits below .01, three digits when two digits would read .05 and
    the value is not .050 itself, >.99 at the top, no leading zero. Carried over unchanged from
    the tables lane's p_jama, which reproduced every P and q cell of Supp Tables 8, 10, 11 and
    14 from their unchanged sources. It was already half up, so nothing here moves it."""
    if v is None or v != v:
        return ""
    if v < 0.001:
        return "<.001"
    s3 = Decimal(f"{float(v):.3f}")
    if s3 < Decimal("0.01"):
        return strip0(str(s3))
    s2 = s3.quantize(Decimal("0.01"), rounding=ROUND_HALF_UP)
    if s2 == Decimal("0.05") and s3 != Decimal("0.050"):
        return strip0(str(s3))
    if s2 >= Decimal("1.00"):
        return ">.99"
    return strip0(str(s2))


def p3(v):
    """Supp Table 9's convention: <.001, else three decimals half up, no leading zero."""
    if v is None or v != v:
        return "NA"
    if v < 0.001:
        return "<.001"
    return strip0(half_up(v, 3))


def pct(v, dp=1):
    return half_up(v, dp)
