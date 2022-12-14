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
            obj.T = [R, p; zeros(1,2), 1];
        end
        
        function obj = assignLink(obj, totalMass, totalHeight, gender, plane)
            switch gender
                case 'F'
                    obj.m = obj.m_ratio_F*totalMass;
                    obj.L = obj.L_ratio_F*totalHeight;
                    obj.c = obj.c_ratio_F*obj.L;
                    switch plane
                        case 'sgt'
                            obj.j = obj.m*(obj.L*obj.r_ratio_sgt_F)^2;
                        case 'frt'
                            obj.j = obj.m*(obj.L*obj.r_ratio_frt_F)^2;
                    end
                case 'M'
                    obj.m = obj.m_ratio_M*totalMass;
                    obj.L = obj.L_ratio_M*totalHeight;
                    obj.c = obj.c_ratio_M*obj.L;
                    switch plane
                        case 'sgt'
                            obj.j = obj.m*(obj.L*obj.r_ratio_sgt_M)^2;
                        case 'frt'
                            obj.j = obj.m*(obj.L*obj.r_ratio_frt_M)^2;
                    end
            end
        end
    end
end
