function hogFeats = computeHog(mag,ang)

ang = wrapToPi(ang);
cellSize = 6;
nBins = 9;
angBins = linspace(-pi/2,pi/2,nBins+1);
hogFeats = zeros(size(mag,1),size(mag,2),size(mag,3));
for u = 1:cellSize:size(mag,2)
    for v = 1:cellSize:size(mag,1)-cellSize
        magCell = mag(v+(1:cellSize),u+(1:cellSize));
        angCell = ang(v+(1:cellSize),u+(1:cellSize));
        
        % compute histogram
        for iBin = 1:nBins
            hogFeats(v,u,iBin) = sum(sum(magCell(angCell>angBins(iBin) & angCell<angBins(iBin+1))));
        end
            
    end
end