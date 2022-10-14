# This file contains the various correlation kernels to use within buildCovariance

function GaspariCohn(r,l)
# Gaspari-Cohn correlation kernel. Takes correlation length (l) and distance between two pts (r) as inputs
if r < l
    b = -(1/4)*(r/l)^5 + (1/2)*(r/l)^4 +(5/8)*(r/l)^3 - (5/3)*(r/l)^2 + 1
elseif r < 2*l
    b = (1/12)*(r/l)^5 - (1/2)*(r/l)^4 +(5/8)*(r/l)^3 + (5/3)*(r/l)^2 - 5*(r/l) + 4 - (2/3)*(l/r)
else
    b = 0
end

end