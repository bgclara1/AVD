function stressCompliant = checkYieldStress(directStress,tensileYieldStress)
    stressCompliant = all(abs(directStress) <= tensileYieldStress);
end



