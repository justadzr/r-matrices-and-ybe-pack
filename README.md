### A small package for classical r-matrices and quantum R-matrices in affine untwisted sl(n)
---
I'm now having trouble to write down the standard part of R_{T, s} given a triple and a continuous datum. Technically speaking, here s can be selected as an element of sl(n). But it troubles me when I tried to produce the standard part. That is, the method I use produces a different standard part in contrast to Polishchuk's formula, although our nonstandard parts are the same (for associative matrices), which is not surprising.

Update: ok I figured out why. I saw the part with e_{-α}\otimes e_{α} as the standard part, and multiplied both sides of it with q^s, which is not what I'm supposed to do. After getting this part out and multiplying the true standard part with q^s, we get the same standard part, and there is no need to think about where s lives.

---
Another update: took the nonassociative triple [3, 4, 0, 0]. All possible passing orders for α1+α2 -> α3+α4 (0, 1/2, 1, 3/2, 2, 5/2, and 3) didn't work in our frame work.

---

TODO: documentations, etc. Everything happens in `test.py`. Some modules contain circular imports.

<del> The current status: I wrote a very naive formula and tested it with some simple nonassociative triples. It seems that the formula worked, but I can't be sure because</del>
1. <del>The passing order was never dealt with rigorously.</del>
2. <del>I did not compute the actual formula of the affine twist.</del>
3. <del>The nonassociative triples I used were way too simple. It might be false if we use triples where T can be applied to a root more than once.</del>
4. <del>I don't know what the "correct" R-matrices should be for these triples, and thus I couldn't verify the formula. </del>

<del>OK so the formula is false for the following triple:
α₁ -> α₄
α₂ -> α₅
α₄ -> α₈
α₅ -> α₇</del>