### A small package for classical r-matrices and quantum R-matrices in affine untwisted sl(n)
---
TODO: documentations, etc. Everything happens in `test.py`. Some modules contain circular imports.
---
<del> The current status: I wrote a very naive formula and tested it with some simple nonassociative triples. It seems that the formula worked, but I can't be sure because</del>
1. <del>The passing order was never dealt with rigorously.</del>
2. <del>I did not compute the actual formula of the affine twist.</del>
3. <del>The nonassociative triples I used were way too simple. It might be false if we use triples where T can be applied to a root more than once.</del>
4. <del>I don't know what the "correct" R-matrices should be for these triples, and thus I couldn't verify the formula. </del>

OK so the formula is false for the following triple:
α₁ -> α₄
α₂ -> α₅
α₄ -> α₈
α₅ -> α₇