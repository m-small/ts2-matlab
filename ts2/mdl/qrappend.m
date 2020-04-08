function [Q,R] = qrappend(Q,R,x)
%  [Q,R] = qrappend(Q,R,x) append a column to a QR factorization.
%	If [Q,R] = qr(A) is the original QR factorization of A,
%	then [Q,R] = qrappend(Q,R,x) changes Q and R to be the
%	factorization of the matrix obtained by appending an extra
%	column, x, to A.
%

[n,m] = size(Q);
if m == 0
  [Q,R] = qr(x,0);
  return;
end

m = m+1;
r = Q'*x; 				% best fit of x by Q
R(:,m) = r; 				% add coeff to R

if 1,% m<n
  q= x - Q*r; 				% q is orthogonal part of x
  f= norm(q);
  R(m,:) = zeros(1,m); 			% update R for q
  R(m,m) = f; 				% f is coeff of q when normalized
else
 % q= x - Q*r; 				% q is orthogonal part of x

  f= norm(q);
 % R(m,:) = zeros(1,m); 			% update R for q
 % R(m,m) = f; 				% f is coeff of q when normalized

end

Q(:,m) = q/f; 				% extend basis by normalized q
