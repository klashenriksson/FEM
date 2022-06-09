function triangle = GetTriangle(edge,t)
    A = find(any(t == edge(1)));
    B = find(any(t == edge(2)));
    triIndex = intersect(A,B);
    triangle = t(:,triIndex);

end