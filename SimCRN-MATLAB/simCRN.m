classdef simCRN < handle
    properties (SetAccess=private)
        ids=cell(1,0);%cell array containing the names of all species, in order of their ID number
        rxns=zeros(0,0);%matrix of reaction rate constants, one column for every reaction, one row for every species, reversible reactions entered in matrix as two separate irreversible reactions
        rxnSpecies=cell(1,0);%cell array containing vectors with the ID's of the reactants in every reaction
        conc=zeros(0,1);%matrix of concentrations for each (species,time)
        time=zeros(1,1);%row vector of timepoints that have been simulated thus far (initially only t=0)
    end
    properties
        lineWidth=2;
        fontSize=21;
    end
    methods
        function this=simCRN()
            %constructor
        end
        function str = printCRN(this)
            for speciesToPrint=1:length(this.ids)
                str='';
                for rxn=1:size(this.rxns,2)
                    rateConstant=this.rxns(speciesToPrint,rxn);
                    str=[str,' ']; %#ok<AGROW>
                    if rateConstant>0 && ~isempty(str)
                        str=[str,'+']; %#ok<AGROW>
                    end
                    if rateConstant~=0
                        str=[str,num2str(rateConstant)]; %#ok<AGROW>
                        for reactantId=this.rxnSpecies{rxn}
                            str=[str,'*',this.ids{reactantId}]; %#ok<AGROW>
                        end
                    end
                end
                str=['d',this.ids{speciesToPrint},'/dt =',str]; %#ok<AGROW>
                disp(str);
            end
            
        end
        function [times,concentrations]=getSpeciesConcentrations(this,speciesName)
            times=this.time;
            concentrations=this.conc(this.getSpeciesIdsByNames(speciesName),:);
        end
        function plotSpecies(this,speciesToPlot,legendOn)%plot all species named in speciesToPlot
            figure
            plot(this.time,this.conc(this.getSpeciesIdsByNames(speciesToPlot),:),'LineWidth',this.lineWidth);
            if legendOn
                legend(speciesToPlot);
            end
            set(gca, 'FontSize',this.fontSize);
        end
        function runSim(this,timeSpan,odeSolverHandle,maxStep,showProgress)%can set maxStep=0 to disable maxStep setting
            for i=1:length(this.rxnSpecies)
            end
            if nargin<3
                odeSolverHandle=@ode15s;%default to the stiff solver
            end
            if nargin<4 || maxStep==0
                opts = odeset();
            else
                opts = odeset('MaxStep',maxStep);
            end
            if nargin<5 || showProgress==1
                opts = odeset(opts,'OutputFcn',@odeprog,'Events',@odeabort);
            end
            [times,concentrations] = odeSolverHandle(@this.getDerivatives,[0,timeSpan],this.conc(:,end),opts);
            this.conc=[this.conc,concentrations'];
            this.time=[this.time,this.time(end)+times'];
        end
        function derivs=getDerivatives(this,t,concentrations) %#ok<INUSL>
            rxnProducts=zeros(1,length(this.rxnSpecies));
            for rxn=1:length(rxnProducts)
                rxnProducts(rxn)=prod(concentrations(this.rxnSpecies{rxn}));
            end
            derivs=zeros(length(concentrations),1);
            for species=1:length(derivs)
                derivs(species)=sum(this.rxns(species,:).*rxnProducts);
            end
        end
        function addRxn(this,reactants,products,kf,kr)
            reactantIds=this.getSpeciesIdsByNames(reactants);
            productIds=this.getSpeciesIdsByNames(products);
            this.rxns=[this.rxns,zeros(size(this.rxns,1),1)];
            this.rxnSpecies{end+1}=reactantIds;
            for i=1:length(reactantIds)
                this.rxns(reactantIds(i),end)=this.rxns(reactantIds(i),end)-kf;
            end
            for i=1:length(productIds)
                this.rxns(productIds(i),end)=this.rxns(productIds(i),end)+kf;
            end
            if kr~=0%add in the reverse reaction
                this.rxns=[this.rxns,zeros(size(this.rxns,1),1)];
                this.rxnSpecies{end+1}=productIds;
                for i=1:length(productIds)
                    this.rxns(productIds(i),end)=this.rxns(productIds(i),end)-kr;
                end
                for i=1:length(reactantIds)
                    this.rxns(reactantIds(i),end)=this.rxns(reactantIds(i),end)+kr;
                end
            end
        end
        function setConcentration(this,varargin)%varargin contains pairs of: 'SpeciesName',concentration. If an odd number of arguments are given, the final argument is taken as the positin in the time vector to override (if even, the time point defaults to the latest timepoint)
            if mod(length(varargin),2)
                timePoint=varargin{end};
            else
                timePoint=size(this.conc,2);
            end
            for i=1:2:length(varargin)
                this.conc(this.getSpeciesIdsByNames(varargin{i}),timePoint)=varargin{i+1};
            end
        end
        function speciesIds=getSpeciesIdsByNames(this,names)
            if iscell(names)
                speciesIds=zeros(1,length(names));
            else
                speciesIds=0;
            end
            for i=1:length(speciesIds)
                if iscell(names)
                    name=names{i};
                else
                    name=names;
                end
                speciesFound=false;
                for id=1:length(this.ids)
                    if strcmp(this.ids{id},name)%species already in list
                        speciesIds(i)=id;
                        speciesFound=true;
                        break
                    end
                end
                if ~speciesFound%species not in list yet...
                    this.ids{length(this.ids)+1}=name;%add it to the list
                    speciesIds(i)=length(this.ids);
                    this.rxns=[this.rxns;zeros(1,size(this.rxns,2))];
                    this.conc(end+1,end)=0;
                end
            end
        end
    end
end