### A small package for classical r-matrices and quantum R-matrices in affine untwisted sl(n)
---
I projected my formula to sl(n) (though some efforts spent on the projection codes...). And it checks out with BD's original formula! It remains to summarize the precise formula of the standard part and put them on Overleaf (I know what the precise formula of the nonstandard part is).

---
Now my program works. I've tested it with all triples of sl(n) for n <= 5. I have some trouble getting the correct standard part out of this formula. Indeed, the nonstandard part of my formula corresponds to the original formula of BD, however the standard parts are vastly different. The same is true for Polishchuk's formula. I suspect my formula is gauge equivalent to BD's.

---
I'm now having trouble to write down the standard part of R_{T, s} given a triple and a continuous datum. Technically speaking, here s can be selected as an element of sl(n). But it troubles me when I tried to produce the standard part. That is, the method I use produces a different standard part in contrast to Polishchuk's formula, although our nonstandard parts are the same (for associative matrices), which is not surprising.

Update: ok I figured out why. I saw the part with e_{-α}\otimes e_{α} as the standard part, and multiplied both sides of it with q^s, which is not what I'm supposed to do. After getting this part out and multiplying the true standard part with q^s, we get the same standard part, and there is no need to think about where s lives.

---
Another update: took the nonassociative triple [3, 4, 0, 0]. All possible passing orders for α1+α2 -> α3+α4 (0, 1/2, 1, 3/2, 2, 5/2, and 3) didn't work in our frame work. In Polishchuk components involving negative roots in \Gamma_1 are treated separately, and he rotates the diagram so that the affine root is not in \Gamma_2. Maybe this is something deserves some thoughts. But thinking about how our formula for [2, 3, 0] coincides with his, it might not be something interesting.

Still, the priority is to understand the affine twist we need. I don't know how to write it down. Might need some help.

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