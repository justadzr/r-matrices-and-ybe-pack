TODO: documentations, etc. Everything happens in `test.py`. Some modules contain circular imports.
---
The current status: I wrote a very naive formula and tested it with some simple nonassociative triples. It seems that the formula worked, but I can't be sure because
1. The passing order was never dealt with rigorously.
2. I did not compute the actual formula of the affine twist.
3. The nonassociative triples I used were way too simple. It might be false if we use triples where T can be applied to a root more than once.
4. I don't know what the "correct" R-matrices should be for these triples, and thus I couldn't verify the formula.