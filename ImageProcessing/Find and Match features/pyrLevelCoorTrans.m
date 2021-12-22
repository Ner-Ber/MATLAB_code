function p_new = pyrLevelCoorTrans(pos, inputlevel, outputlevel)
%pyrLevelCoorTrans transfoms coordinates from one level of te pyramid to
%another

    p_new = (2^(inputlevel - outputlevel))*(pos - 1) +1;

end