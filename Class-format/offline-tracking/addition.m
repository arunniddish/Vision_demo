% classdef untitled
%     %UNTITLED Summary of this class goes here
%     %   Detailed explanation goes here
% 
%     properties
%         Property1
%     end
% 
%     methods
%         function obj = untitled(inputArg1,inputArg2)
%             %UNTITLED Construct an instance of this class
%             %   Detailed explanation goes here
%             obj.Property1 = inputArg1 + inputArg2;
%         end
% 
%         function outputArg = method1(obj,inputArg)
%             %METHOD1 Summary of this method goes here
%             %   Detailed explanation goes here
%             outputArg = obj.Property1 + inputArg;
%         end
%     end
% end

classdef addition

    properties
        a;
        b;
        c;
    end 

    methods
        function obj = addition(paramet)
            for i = 1:5
            obj.a = paramet.aa;
            obj.b = paramet.bb;
            h(i) = summ(obj);
            end
            obj.c = h;
        end

        function output = summ(obj)
            output = obj.a + obj.b;
        end
    end
end

