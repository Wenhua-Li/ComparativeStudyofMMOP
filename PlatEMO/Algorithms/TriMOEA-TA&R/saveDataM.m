function saveDataM(obj)
                folder = fullfile('Data',func2str(obj.algorithm));
                [~,~]  = mkdir(folder);
                Population     = obj.result{end};
                PF             = obj.PF;
                metric.runtime = obj.runtime;
                name = strcat('%s_%s_M%d_%d');
                addre= strcat(name,'.mat');
                save(fullfile(folder,sprintf(addre,func2str(obj.algorithm),func2str(obj.problem),obj.M,obj.run)),'Population','PF','metric');
end