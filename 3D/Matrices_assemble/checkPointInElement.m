function isInElement = checkPointInElement(xi, range)
% Check if the point is in the element
isInElement = (range(2) >= xi(1) && xi(1) >= range(1)) && (range(4) >= xi(2) && xi(2) >= range(3));
end