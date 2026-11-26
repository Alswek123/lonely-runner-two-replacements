# Four replacements of consecutive speeds

This directory preserves the earlier proof, source and replay data. It classifies four arbitrary replacements for all total runner counts $n\ge9$, and records the known classifications at $n=5,6$. Its original cutoff on the smallest new speed at $n=7,8$ was removed in the subsequent equality paper. Together the two drafts cover every size. [Subsequent proof, Corollary 2](../equality/equality.pdf) · [Earlier proof](proof.pdf) · [Original Chinese notes and full rebuild instructions](README.zh-CN.md)

For $n=9$ the only result is $2\{1,\ldots,8\}$. For $n\ge10$ the draft characterizes the tight results using the simultaneous Goddyn–Wong gcd conditions. The first such size is $148431434$. A companion arithmetic converse supplies the implication from the geometric pairing to the gcd conditions. [Arithmetic converse](paired-converse.pdf)

Run `python verify.py` for the arithmetic checks. For full regeneration, including the interval covering census, run `python verify.py --full` with Python 3.9+ and a C++17 compiler with 128-bit integer support. A multiword compiler command is accepted, for example `python verify.py --full --compiler zig c++`. Run Python with its default assertions enabled.
