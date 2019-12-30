function I=I_note(x,inf_)
k=[1:inf_]';
I=1+sum((((x./2).^k)./factorial(k)).^2);
end