function [g] = GlottalConstraint(N,repi,iepi)
%repi - Epilaryngeal weight for range of vocal tract region (Story et.al
%used 6)
%iepi - Epilayngeal sections above glottis (Story et.al used 5)

for n = 1:N
    if (n>=1 && n<=iepi-1)
        g(n) = 0;
    elseif (n>iepi && n<=N)
        g(n) = 1 - exp(-log(16).*((n - iepi)./repi).^2);
    end
end
g = g.';
end