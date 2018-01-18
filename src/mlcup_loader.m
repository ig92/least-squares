function Xtilde = mlcup_loader(path, varargin)
    Xtilde = csvread(path);
    Xtilde = Xtilde(1:end, 2:11);
    Nargs = size(varargin,2);
    if (Nargs ~= 0)
       for i = 1 : Nargs
           f = varargin(i);
           Xtilde = f(Xtilde);
       end
    end
end