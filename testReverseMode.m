function testReverseMode
    Extra.g = @metricSphere;
    
    AD_fun = ADfun(@testFun,4);
    options = setopt('revprod',[1 0;0 1;0 0;0 0]);
    x = [pi/3;pi/3];
    [c, c_i] = feval(AD_fun, x, Extra, options)
    g_i = metricSphere_i(x)
end

function z = testFun(x,Extra)
    z = reshape(x,3,[]);
end