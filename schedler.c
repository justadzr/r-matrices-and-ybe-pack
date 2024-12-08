/*  ggstest2.c  GGS test, version 2
    Copyright 1999 by Travis Schedler, <schedler@fas.harvard.edu>.  
    All rights reserved.
    This program is licensed under the GNU General Public License.
    Verbatim distribution is encouraged.  This program has no warranty.


    Instructions on usage: Create a text file called "triples" in the
    directory with this program.  Once the file is created, the program 
    should be compiled and run with no arguments, with output given in
    standard output.  The format of the output is described in the
    output section.  The file "triples" should be formatted as follows:

    Note: the commas and periods below are essential!  a_i, b_i, etc. 
    are to be filled in by natural numbers.

    Line 1:  a_1,a_2,...,a_k.
    Line 2:  b_1,b_2,...,b_l.
            .
            .
            .

    Each line contains one triple.  In this case, the first triple is in 
    sl(k-1), the second sl(l-1).  We must have k,l < n.  The notation is 
    this: a_i indicates \tau(\alpha_i) = \alpha_{a_i}, unless a_i = 0, 
    in which case \alpha_i is not in \Gamma_1.

    The commas and periods are essential and must be as indicated!
    That is, commas separate n,d and elements of a given triple, and
    periods separate n,d and the triples from each other.  The line
    breaks need not be as indicated, however.
  


    Example "triples" file:

    2,3,4,0.
    3,4,0,1.

    This file indicates n <= 5, and 100 are enough terms in all
    polynomials.  The first line is the Cremmer-Gervais triple,
    \alpha_1 -> \alpha_2, \alpha_2 -> \alpha_3, \alpha_3 ->
    \alpha_4. So \alpha_4 is not in \Gamma_1 in this case.  The second
    line is the "generalized" Cremmer-Gervais triple, \alpha_1 ->
    \alpha_3, \alpha_2 -> \alpha_4, \alpha_4 -> \alpha_1.  Here,
    \alpha_3 is not in \Gamma_1.

    
    Instructions on interpreting the output: The output is given to
    the standard output.  First, the program describes the triple it
    has been fed, in a manner such as 1 -> 2, 2 -> 3, which describes
    \tau(\alpha_1) = \alpha_2; \tau(\alpha_2) = \alpha_3.  Then, the
    program prints the \tilde r^0 and R matrices.  These are printed
    in a sequence of k matrices, each k x k, for the space sl(k).  The
    i-th matrix is indexed by j and k, and prints the entries a_{ik}^j
    e_{ij} \otimes e_{k,i+k-j}.  Thus, the k-th entry on the j-th line
    of the i-th matrix is a_{ik}^j e_{ij} \otimes e_{k,i+k-j}.  For
    \tilde r^0, the numbers should all be divided by the given
    denominator, and for R, the exponents of q should all be divided
    by the given denominator.


    Notes on compilation: I recommend using gcc; I always have good
    experience with the -O2 option.  I have only tested this program
    with gcc, on x86, UltraSparc, and R4600 environments running
    Linux 2.0.x, SunOS 5.6, and IRIX 6.2, respectively.


    If, for some reason, memory is a problem for you and you would like
    a version that does not store all of R in memory at one time, let me
    know.  Send any other comments or requests to me as well.
    Good luck with the program!  
*/

#include <stdio.h>
#include <stdlib.h>

#define dm 20
#define dmp 500

struct poly {
int posdeg;
int negdeg;
int pos[dmp];
int neg[dmp];
int con;
};

struct sp {
int c;
int a;
int b;
};

int rot[dm][dm], a[dm][dm][dm], c[dm][dm][dm], eps[dm][dm][dm];
struct sp r[dm][dm][dm];


int gcd(int a, int b) {
  int temp = abs(a), temp2 = abs(b), done;
  if(temp == 1 || temp2 == 1) done = 1;
  else if(temp == 0 || temp2 == 0) done = temp+temp2;
  while(temp > 1 && temp2 > 1) {
    temp = temp % temp2;
    if(temp == 1) done = 1;
    else if(temp == 0) done = temp2;
    else {
      temp2 = temp2 % temp;
      if(temp2 == 1) done = 1;
      else if(temp2 == 0) done = temp;
    }
  }
  return(done);
}

int gcdm(int mat[dm][dm], 
		      int new[dm][dm], int size, 
		      int factor) {
  int temp = abs(mat[0][0]), ret, count, count2;
  for(count = 0; count < size && temp != 1 ; count++) {
    for(count2 = 0; count2 < size-1 && temp != 1; count2++)
      temp = gcd(temp, mat[count][count2+1]);
    if(count < size-1) temp = gcd(temp, mat[count+1][0]);
  }
  temp = gcd(temp, factor);
  ret = factor/temp;
  for(count = 0; count < size; count++) {
    for(count2 = 0; count2 < size; count2++) {
      new[count][count2] = mat[count][count2]/temp;
    }
  }
  return(ret);
}

void setupt(int tau[dm], int iti[dm][dm], 
	    int sg) {
  int count, count2;
  for(count = 0; count < sg; count++) {
    for(count2 = 0; count2 < sg; count2++) {
      iti[count][count2]=0;
    }
  }
  for(count = 0; count < sg; count++) {
    count2 = count;
    while(count2 != sg) {
      iti[count2][count]++;
      count2 = tau[count2];
    }
  }
}

int adp(int k, int l, 
		     int sg) {
  int ans;
  if(k == sg || l == sg) ans = 0;
  else if(abs(k-l)==1) ans = -1;
  else if(k==l) ans = 2;
  else ans = 0;
  return(ans);
}

void setupb(int tau[dm], int b[dm][dm], 
	    int sg) {
  int count, count2;
  for(count = 0; count < sg; count++) {
    b[count][count]=0;
    if(tau[count] != sg) {
      for(count2 = 0; count2 < sg; count2++) {
	b[count][count2]=adp(count,count2,sg)+adp(tau[count],count2,sg)-
	  adp(count,tau[count2],sg)-adp(tau[count],tau[count2],sg);
	b[count2][count]=0-b[count][count2];
      }
    } else 
      for(count2 = 0; count2 < count; count2++)
	if(tau[count2] == sg) {b[count][count2] = 0; b[count2][count] = 0;}
    
  }
}

void altoe(int al[dm][dm], int e[dm][dm], 
	   int sg) {
  int count, count2, temp[dm][dm];
  for(count = 0; count < sg; count++) {
    temp[count][0]=0;
    for(count2 = 0; count2 < sg; count2++) {
      temp[count][0]+=(sg-count2)*al[count][count2];
    }
    for(count2 = 1; count2 < sg+1; count2++) {
      temp[count][count2]=-(sg+1)*al[count][count2-1]+temp[count][count2-1];
    }
  }
  for(count = 0; count < sg+1; count++) {
    e[0][count]=0;
    for(count2 = 0; count2 < sg; count2++) {
      e[0][count]+=(sg-count2)*temp[count2][count];
    }
    for(count2 = 1; count2 < sg+1; count2++) {
      e[count2][count]=-(sg+1)*temp[count2-1][count]+e[count2-1][count];
    }
  }
}
    

void trans(int a[dm][dm], int b[dm][dm], 
	   int size) {
  int count, count2;
  for(count = 0; count < size; count++) {
    b[count][count] = a[count][count];
    for(count2 = 0; count2 < count; count2++) {
      b[count][count2] = a[count2][count];
      b[count2][count] = a[count][count2];
    }
  }
}

void mult(int a[dm][dm], int b[dm][dm], 
	  int c[dm][dm], int size) {
  int count, count2, count3, temp;
  for(count = 0; count < size; count++) {
    for(count2 = 0; count2 < size; count2++) {
      temp = 0;
      for(count3 = 0; count3 < size; count3++) {
	temp += a[count][count3]*b[count3][count2];
      }
      c[count][count2] = temp;
    }
  }
}

void printmat(int mat[dm][dm], int size) {
  int count, count2, tempn;
  for(count = 0; count < size; count++) {
    for(count2 = 0; count2 < size; count2++) {
      if(mat[count][count2] >= 0) printf(" ");
      tempn = abs(mat[count][count2]);
      if(!tempn) tempn=1;
      while(tempn < 10000) {
	tempn *= 10;
	printf(" ",tempn);
      }
      printf("%d\t",mat[count][count2]);
    }
    printf("\n");
  }
  printf("\n");
}

void linear(int tau[dm], int sg, 
	    int *d) { 
  int open = 1, count, b[dm][dm], iti[dm][dm], 
    count2, itit[dm][dm], temp1[dm][dm], temp2[dm][dm], tempn;
  setupt(tau, iti, sg);
  trans(iti, itit, sg);
  setupb(tau, b, sg);
  mult(itit, b, temp1, sg);
  mult(temp1, iti, temp2, sg);
  altoe(temp2, temp1, sg);
  *d = gcdm(temp1, rot, sg+1, 2*(sg+1)*(sg+1));
  printf("\nrtilde, fractions of denominator %d:\n",*d);
  printmat(rot, sg+1);
}

int red(int x, int n) {
  int t = x;
  while(t < 0) t+=n;
  while(t > n-1) t-=n;
  return t;
}


void spprint(struct sp *a) {
  if((a->c) == 0) printf("0"); 
  else {
    if(abs(a->c) == 1) {
      if(a->c == -1) printf("-");
      if((a->a) == (a->b) && (a->a) == 0) printf("1");
    }
    else printf("%d", a->c);
    if(abs(a->c) > 1 && ((a->a) || (a->b))) printf(" ");
    if((a->a) != 0) printf("q^%d", a->a);
    if((a->a) && (a->b)) printf(" ");
    if((a->b) != 0) printf("qh^%d", a->b);
  }
}

void reduce(struct poly *a) {
  int count;
  for(count = a->posdeg-1; count >= 0 && a->pos[count] == 0; count--);
  a->posdeg = count+1;
  for(count = a->negdeg-1; count >= 0 && a->neg[count] == 0; count--);
  a->negdeg = count+1;
}

void polqpadd(struct poly *a, int c, int p) {
  int count;
  if(p == 0) a->con += c;
  else if(p > 0) {
    if(a->posdeg < p) {
      for(count = a->posdeg; count < p - 1; count++)
	a->pos[count] = 0;
      a->posdeg = p;
      a->pos[p - 1] = c;
    } else a->pos[p - 1] += c;
  } else {
    if(a->negdeg < (-p)) {
      for(count = a->negdeg; count < -1 - p; count++)
	a->neg[count] = 0;
      a->negdeg = (-p);
      a->neg[-1 - p] = c;
    } else a->neg[-1 - p] += c;
  }
  reduce(a);
}

void sppoladd(struct sp *a, struct poly *b, int d) {
  int count;
  if(a->c != 0) {
    if(a->b == 0) { 
      polqpadd(b, a->c, a->a);
    } else if(a->b == 1) {
      polqpadd(b, a->c, d+(a->a));
      polqpadd(b, 0-(a->c), (a->a)-d); 
    } else if(a->b == 2) {
      polqpadd(b, a->c, (d << 1)+(a->a));
      polqpadd(b, a->c, (a->a)-(d << 1)); 
      polqpadd(b, -((a->c) << 1), a->a);
    } else if(a->b == 3) {
      polqpadd(b, a->c, (3*d)+(a->a));
      polqpadd(b, -(a->c), (a->a)-(3*d)); 
      polqpadd(b, 3*(a->c), (a->a)-d);
      polqpadd(b, (-3)*(a->c), d+(a->a)); 
    }
  }
}

void zero(struct poly *a) {
  a->posdeg = 0;
  a->negdeg = 0;
  a->con = 0;
}

void spscalmult(struct sp *a, int mult) {
  if(mult == 0) {a->b = 0; a->a = 0; a->c = 0;}
  else a->c *= mult;
}

void spass(struct sp *a, struct sp *b) {
  b->a = a->a;
  b->b = a->b;
  b->c = a->c;
}

void spmult(struct sp *a, struct sp *b, struct sp *c) {
  c->c = b->c * a->c;
  c->a = b->a + a->a;
  c->b = b->b + a->b;
}

int readline(FILE *silly, int tau[dm], 
			  int *n) {
  int count, count2, good = 1, ch, temp[dm], tempn, sign;
  *n = 1;
  while(good == 1) {
    tau[*n-1] = 0;
    while(good == 1) {
      ch = fgetc(silly);
      if(ch == EOF) good = 0;
      else if(ch == ',') good = 3;
      else if(ch == '.') good = 2;      
      else if(ch>47 && ch<58) {
	tau[*n-1] *= 10;
	tau[*n-1] += ch - 48;
      }
    }
    if(good == 3) {
      tau[*n-1]--;
      (*n)++;
      good = 1;
    }
  }
  tau[*n-1]--;
  (*n)++;
  for(count = 0; count < *n - 1; count++) 
    if(tau[count] < 0) tau[count] = (*n) - 1;
  return(good);
}


void tens2mult(int a[dm][dm][dm], int b[dm][dm][dm], 
	       int c[dm][dm][dm], int size) {
  int count, count2, count3, count4;
  for(count = 0; count < size; count++) {
    for(count2 = 0; count2 < size; count2++) {
      for(count3 = 0; count3 < size; count3++) {
        c[count][count2][count3] = 0;
	for(count4 = 0; count4 < size; count4++) {
	  if(count+count3-count4 == red(count+count3-count4,size) &&
	     count+count3-count2 == red(count+count3-count2,size)) {
	    c[count][count2][count3] += a[count][count4][count3] * 
	      b[count4][count2][count+count3-count4]; 
	  }
	}
      }
    }
  }
}

void tens2add(int a[dm][dm][dm], int b[dm][dm][dm], int size) {
  int count, count2, count3;
  for(count = 0; count < size; count++) {
    for(count2 = 0; count2 < size; count2++) {
      for(count3 = 0; count3 < size; count3++) {
        b[count][count2][count3] += a[count][count2][count3];
      }
    }
  }
}

void scmult(int k, int a[dm][dm][dm], int size) {
  int count, count2, count3;
  for(count = 0; count < size; count++) {
    for(count2 = 0; count2 < size; count2++) {
      for(count3 = 0; count3 < size; count3++) {
        a[count][count2][count3] *= k;
      }
    }
  }
}

void setupa(int tau[dm], int n) {
  int count, count2, count3, temp[dm], s, f;
  for(count = 0; count < n; count++)
    for(count2 = 0; count2 < n; count2++)
      for(count3 = 0; count3 < n; count3++)
	  a[count][count2][count3]=0;
  
  for(count = 0; count < n-1; count++) {
    for(count2 = count; count2 < n-1; count2++) 
      temp[count2]=tau[count2];
    while(temp[count] < n-1) {
      for(count2 = count; count2 < n-1 && temp[count2] < n-1; count2++) {
	if(temp[count2] > temp[count]) {
	  s = temp[count]; f = temp[count2];
	  a[count2+1][count][s]++;
	  a[s][f+1][count2+1]--;
	}
	else {
	  s = temp[count2]; f = temp[count];
          if((count2 - count) % 2 == 0) {
            a[count2+1][count][s]++;
            a[s][f+1][count2+1]--;
          } else {
            a[count2+1][count][s]--;
            a[s][f+1][count2+1]++;
          }
	}
      }
      for(count2 = count; count2 < n-1 && temp[count2] < n-1; count2++)
	temp[count2]=tau[temp[count2]];      
    }
  }
}

void setupc(int n) {
  int count, count2, count3;
  for(count = 0; count < n; count++)
    for(count2 = 0; count2 < n; count2++)
      for(count3 = 0; count3 < n; count3++)
	  c[count][count2][count3]=0;
  
  for(count = 0; count < n; count++) {
    for(count2 = 0; count2 < count; count2++) {
      c[count2][count][count]--;
      c[count][count2][count2]++;
    }
  }
}

void retr(struct sp *pol, int n, int d, int i, int j, int k) {
  if(i == j) {
    pol->c = 1;
    pol->b = 0;
    if(i == k)
      pol->a = d << 1;
    else
      pol->a = rot[i][k] << 1;
  }
  else {
    pol->b = 1;
    if(j == k) {
      pol->a = 0;      
      if(i > k || i+k<j || i+k>=j+n) pol->c = 1;
      else pol->c = 0;
    } else { 
      pol->a = d*a[i][j][k]*eps[i][j][k]+rot[i][k]+rot[j][i+k-j];
      pol->c = a[i][j][k];
    }
  }
}

void setupr(int n, int d) {
  int ct1, ct2, ct3;
  printf("R, powers have denominator %d:\n",d*2);
  for(ct1 = 0; ct1 < n; ct1++) {
    for(ct2 = 0; ct2 < n; ct2++) {
      for(ct3 = 0; ct3 < n; ct3++) {
	retr(&r[ct1][ct2][ct3], n, d, ct1, ct2, ct3);
	spprint(&r[ct1][ct2][ct3]);
	if(ct3 < n-1) printf(",   ");
      }
      printf("\n");
    }
    printf("\n");
  }
}

int qybe(int n, int d, struct poly *t1, struct sp *t3, struct sp *t4) {
  int i,j,k,l,m,o,good = 1;
  int counter = 0;
  for(i = 0; i < n && good; i++) 
    for(j = 0; j < n && good; j++) 
      for(k = 0; k < n && good; k++) 
        for(l = 0; l < n && good; l++) 
          for(m = 0; m < n && good; m++) { if(m+k+i-j-l == red(m+k+i-j-l,n)) {
            zero(t1);
            for(o = 0; o < n; o++) { 
	      if(k+i-o == red(k+i-o,n) && m+k+i-o-j == red(m+k+i-o-j,n)) {
	      spmult(&r[i][k+i-o][k], &r[k+i-o][j][m], t3);
	      spmult(t3, &r[o][l][m+k+i-o-j], t4);
        counter+=1;
	      sppoladd(t4, t1, d << 1);
	      }

	      if(j+l-o == red(j+l-o,n) && m+k-o == red(m+k-o,n)) {
	      spmult(&r[k][o][m], &r[i][j+l-o][m+k-o], t3);
	      spmult(t3, &r[j+l-o][j][o], t4);
	      spscalmult(t4, -1);
        counter+=1;
	      sppoladd(t4, t1, d << 1);
	      }
	    }
            if(t1->con != 0 || t1->posdeg != 0 || t1->negdeg != 0) good = 0; 
          }
	  }
  printf("%d multiplications are performed", counter);
  return good;
}

int hecke(int n, int d, struct poly *t1, struct sp *t3) {
  int i,j,k,l,m,o,good = 1;
  for(i = 0; i < n && good; i++) 
    for(j = 0; j < n && good; j++) 
      for(k = 0; k < n && good; k++) { if(k+i-j == red(k+i-j,n)) { 
	zero(t1);
	if(i == j) t1->con = -1;
	spass(&r[k][j][i], t3);
	t3->b++;
	spscalmult(t3, -1);
	sppoladd(t3, t1, d << 1);
        for(l = 0; l < n; l++) { if(k+i-l == red(k+i-l,n)) {
	  spmult(&r[k][l][i], &r[k+i-l][j][l], t3);
	  sppoladd(t3, t1, d << 1);
	}
	}
	if(t1->con != 0 || t1->posdeg != 0 || t1->negdeg != 0) good = 0;
      }
      }
  return good;
}

int main() {
  int open = 1, count, count2, tau[dm], n=0, m=0, d=0, 
    temp[dm][dm][dm], temp2[dm][dm][dm], tempn;
  struct poly tm, tm2;
  struct sp tsp, tsp2;
  FILE *silly;

  if ((silly =
       fopen("C:\\Users\\danny\\Desktop\\Temp\\explicit-quantization\\codes\\triples.txt", "r")) == NULL) 
    fprintf(stderr, "Cannot open %s\n", "triples");  

  while(open) {
    m = 0;
    if(open = readline(silly, tau, &n)) {
      for(count = 0; count < n-1; count++) {
	if(tau[count] < n-1) {m++; printf("%d -> %d,",count+1, tau[count]+1);}
      }
      printf("m=%d,n=%d::",m,n);
      linear(tau, n-1, &d);

      tempn = gcd(2,d); 
      d /= tempn;
      tempn = 2/tempn;
      for(count = 0; count < n; count++) {
	for(count2 = 0; count2 < n; count2++) {
	  rot[count][count2] *= tempn;
	}
      }

      setupa(tau, n);
      setupc(n);

      tens2mult(a, c, temp, n);
      tens2mult(c, a, temp2, n);
      tens2add(temp2, temp, n);
      tens2mult(a, a, eps, n);
      scmult(2, eps, n);
      tens2add(temp, eps, n);

      setupr(n, d); 
      if(hecke(n, d, &tm, &tsp)) 
	printf("(hecke passed)"); else printf("(hecke failed)");
      if(qybe(n, d, &tm, &tsp, &tsp2)) 
      printf("(qybe passed)"); else printf("(qybe failed)"); 
      printf("\n\n");
      printf("%d", n);
    } 
  }
  fclose(silly);   
}