function x = parametric_curve(s,type)

switch type

    case 'bubble'

        [x,y] = bubble(s);

    case 'ellipse'

        [x,y] = ellipse(s,0,a,b);

    otherwise
        error('not implemented');

end

x = [x y];
