classdef ClassSE2
    properties
        R
        p
        T
    end
    methods
        function obj = ClassSE2(rotationAngle, longitudinalTranslation, transverseTranslation)
            if nargin > 2
                obj.R = rot(rotationAngle);
                obj.p = [longitudinalTranslation; transverseTranslation];
            else
                obj.R = rot(rotationAngle);
                obj.p = [longitudinalTranslation; 0];
            end
            obj.T = [obj.R, obj.p; zeros(1,2), 1];
        end
        
        function obj = mtimes(o1, o2)
            obj.T = o1.T*o2.T;
            obj.p = simplify(obj.T(1:2,3), 'steps',3); 
            obj.R = obj.T(1:2, 1:2);
        end
    end
end
