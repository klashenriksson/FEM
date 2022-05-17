function [An, bn, remainingdofs] = LockDofs(A, b, dofs, values)
    n_total = length(b);
    remainingdofs = (1:n_total);
    remainingdofs(dofs) = [];

    An = A(remainingdofs, remainingdofs);
    bn = b(remainingdofs) - A(remainingdofs, dofs)*values;
end