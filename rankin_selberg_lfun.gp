\\ ============================================================
\\ Tensor product L-function of two elliptic curves E1, E2
\\ Computes:
\\ - Euler factors
\\ - L(E1 x E2, 2) via Dokchitser
\\ ============================================================
\\ *Define elliptic curves (minimal models)
E1 = ellinit([0, 1, 1, -89, 295]);
E2 = ellinit([0, 0, 1, -199, 1092]);

Cond1 = ellglobalred(E1)[1];
Cond2 = ellglobalred(E2)[1];

if(vecmax(factor(Cond1)[,2]) > 1, print("E1 is not semistable"); exit());
if(vecmax(factor(Cond2)[,2]) > 1, print("E2 is not semistable"); exit());

Cond = 1;
facs = factor(Cond1*Cond2)[,1];
for(i = 1, #facs, Cond = Cond * facs[i]^2);

print("Cond1 = ", Cond1);
print("Cond2 = ", Cond2);
print("Tensor conductor = ", Cond);

\\ ------------------------------------------------------------
\\ Euler factors for tensor product at good primes
\\ If ap, bp are traces:
\\ L_p(T) = prod_{i,j}(1 - alpha_i * beta_j T)
\\ = 1 - ap*bp*T + (p*(ap^2 + bp^2 - 2p))*T^2 - ap*bp*p^2*T^3 + p^4*T^4
\\ factors at bad primes are only correct if E1 and E2 are semistable.

tensor_euler_factor(p, a, b) = {
  local(c1,c2,c3,c4);
  c1 = -a*b;
  c2 = p*(a^2 + b^2 - 2*p);
  c3 = -a*b*p^2;
  c4 = p^4;
  1/(1 + c1*x + c2*x^2 + c3*x^3 + c4*x^4)
};

tensor_euler_factor_st_good(p, a, b) = {
  1/(1 - a*b * x + p*x^2)
};

tensor_euler_factor_st_st(p, a, b) = {
  1/(1 - a*b*x)/(1 - p*a*b*x)
};

euler(p,d) = {
  if(Cond2 % p != 0,
    if(Cond1 % p != 0,
      tensor_euler_factor(p, ellap(E1, p), ellap(E2, p)),
      tensor_euler_factor_st_good(p, ellap(E1, p), ellap(E2, p))),
    if(Cond1 % p != 0,
      tensor_euler_factor_st_good(p, ellap(E1,p), ellap(E2,p)),
      tensor_euler_factor_st_st(p, ellap(E1, p), ellap(E2, p))))
};

L = lfuncreate([euler, 0, [-1,0,0,1], 3, Cond, 1]);
print("lfuncheck: ", lfuncheckfeq(L));
print("L'(1) = ", lfun(L, 1, 1));
print("L(2) = ", lfun(L,2));
