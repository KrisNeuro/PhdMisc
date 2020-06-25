function [DeltaPowerCorr] = PlotDeltaRegression(outopt,Rsquare,stage,idx, ...
            il_delta_overtime,il_delta_overtime_1,il_delta_overtime_2,il_delta_overtime_3, ...
            dhip_delta_overtime,dhip_delta_overtime_1,dhip_delta_overtime_2,dhip_delta_overtime_3, ...
            vhip_delta_overtime,vhip_delta_overtime_1,vhip_delta_overtime_2,vhip_delta_overtime_3, ...
            pl_delta_overtime,pl_delta_overtime_1,pl_delta_overtime_2,pl_delta_overtime_3)

% Plot delta power regression
DeltaPowerCorr = figure('Color',[1 1 1],'units','normalized','outerposition',[0 0 1 1]);
        subplot(2,3,1) %IL-VHIP
        if outopt
            plot(il_delta_overtime,vhip_delta_overtime,'+r')
            hold on
        end
        plot(il_delta_overtime_2,vhip_delta_overtime_1,'+k')
        title('IL-vHPC delta power correlation')
        xlabel('IL delta power (mV^2/Hz)')
        ylabel('vHPC delta power (mV^2/Hz)')
        box off

        subplot(2,3,2) %IL-DHIP
        if outopt
            plot(il_delta_overtime,dhip_delta_overtime,'+r')
            hold on
        end
        plot(il_delta_overtime_1,dhip_delta_overtime_1,'+k')
        title('IL-dHPC delta power correlation ')
        xlabel('IL delta power (mV^2/Hz)')
        ylabel('dHPC delta power (mV^2/Hz)')
        box off

        subplot(2,3,3) %DHIP-VHIP
        if outopt
            plot(dhip_delta_overtime,vhip_delta_overtime,'+r')
            hold on
        end
        plot(dhip_delta_overtime_3,vhip_delta_overtime_3,'+k')
        title('dHPC-vHPC delta power correlation')
        xlabel('dHPC delta power (mV^2/Hz)')
        ylabel('vHPC delta power (mV^2/Hz)')
        box off

        subplot(2,3,4) %PL-VHIP
        if outopt
            plot(pl_delta_overtime,vhip_delta_overtime,'+r')
            hold on
        end
        plot(pl_delta_overtime_2,vhip_delta_overtime_2,'k+');
        title('PL-vHPC delta power correlation')
        xlabel('PL delta power (mV^2/Hz)')
        ylabel('vHPC delta power (mV^2/Hz)')
        box off

        subplot(2,3,5) %PL-DHIP
        if outopt
            plot(pl_delta_overtime,dhip_delta_overtime,'+r')
            hold on
        end
        plot(pl_delta_overtime_1,dhip_delta_overtime_2,'k+');
        title('PL-dHPC delta power correlation')
        xlabel('PL delta power (mV^2/Hz)')
        ylabel('dHPC delta power (mV^2/Hz)')
        box off

        subplot(2,3,6) %IL-PL
        if outopt
            plot(il_delta_overtime,pl_delta_overtime,'+r')
            hold on
        end
        plot(il_delta_overtime_3,pl_delta_overtime_3,'k+');
        title('IL-PL delta power correlation')
        xlabel('IL delta power (mV^2/Hz)')
        ylabel('PL delta power (mV^2/Hz)')
        box off
        
        % Create textboxes displaying Rsquared values (outliers removed)
        %IL-VHIP delta
        annotation(gcf,'textbox','String',['r^2= ', num2str(Rsquare.Delta.(stage).il_vhip(idx))],'FitHeightToText','off',...
            'LineStyle','none','Position',[0.2711 0.5888 0.0650 0.0345]);

        % IL-DHIP  delta
        annotation(gcf,'textbox','String',['r^2= ', num2str(Rsquare.Delta.(stage).il_dhip(idx))],'FitHeightToText','off',...
            'LineStyle','none','Position',[0.5538 0.5900 0.0629 0.0345]);
        % DHIP-VHIP delta
        annotation(gcf,'textbox','String',['r^2= ', num2str(Rsquare.Delta.(stage).dhip_vhip(idx))],'FitHeightToText','off',...
            'LineStyle','none',...
            'Position', [0.8441 0.5875 0.0684 0.0319]);
        % PL-VHIP delta
        annotation(gcf,'textbox','String',['r^2= ', num2str(Rsquare.Delta.(stage).pl_vhip(idx))],'FitHeightToText','off',...
            'LineStyle','none',...
            'Position',[0.2649 0.1149 0.0712 0.0370]);
        % PL-DHIP delta
        annotation(gcf,'textbox','String',['r^2= ', num2str(Rsquare.Delta.(stage).pl_dhip(idx))],'FitHeightToText','off',...
            'LineStyle','none',...
            'Position',[0.5448 0.1175 0.0698 0.0307]);
        % IL-PL delta
        annotation(gcf,'textbox','String',['r^2= ', num2str(Rsquare.Delta.(stage).il_pl(idx))],'FitHeightToText','off',...
            'LineStyle','none',...
            'Position',[0.8232 0.1124 0.0747 0.0409]);   
end %function