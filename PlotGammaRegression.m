function [GammaPowerCorr] = PlotGammaRegression(outopt,Rsquare,stage,idx, ...
            il_gamma_overtime,il_gamma_overtime_1,il_gamma_overtime_2,il_gamma_overtime_3, ...
            dhip_gamma_overtime,dhip_gamma_overtime_1,dhip_gamma_overtime_2,dhip_gamma_overtime_3, ...
            vhip_gamma_overtime,vhip_gamma_overtime_1,vhip_gamma_overtime_2,vhip_gamma_overtime_3, ...
            pl_gamma_overtime,pl_gamma_overtime_1,pl_gamma_overtime_2,pl_gamma_overtime_3)

% Plot gamma power regression
GammaPowerCorr = figure('Color',[1 1 1],'units','normalized','outerposition',[0 0 1 1]);
        subplot(2,3,1) %IL-VHIP
        if outopt
            plot(il_gamma_overtime,vhip_gamma_overtime,'+r')
            hold on
        end
        plot(il_gamma_overtime_2,vhip_gamma_overtime_1,'+k')
        title('IL-vHPC gamma power correlation')
        xlabel('IL gamma power (mV^2/Hz)')
        ylabel('vHPC gamma power (mV^2/Hz)')
        box off

        subplot(2,3,2) %IL-DHIP
        if outopt
            plot(il_gamma_overtime,dhip_gamma_overtime,'+r')
            hold on
        end
        plot(il_gamma_overtime_1,dhip_gamma_overtime_1,'+k')
        title('IL-dHPC gamma power correlation ')
        xlabel('IL gamma power (mV^2/Hz)')
        ylabel('dHPC gamma power (mV^2/Hz)')
        box off

        subplot(2,3,3) %DHIP-VHIP
        if outopt
            plot(dhip_gamma_overtime,vhip_gamma_overtime,'+r')
            hold on
        end
        plot(dhip_gamma_overtime_3,vhip_gamma_overtime_3,'+k')
        title('dHPC-vHPC gamma power correlation')
        xlabel('dHPC gamma power (mV^2/Hz)')
        ylabel('vHPC gamma power (mV^2/Hz)')
        box off

        subplot(2,3,4) %PL-VHIP
        if outopt
            plot(pl_gamma_overtime,vhip_gamma_overtime,'+r')
            hold on
        end
        plot(pl_gamma_overtime_2,vhip_gamma_overtime_2,'k+');
        title('PL-vHPC gamma power correlation')
        xlabel('PL gamma power (mV^2/Hz)')
        ylabel('vHPC gamma power (mV^2/Hz)')
        box off

        subplot(2,3,5) %PL-DHIP
        if outopt
            plot(pl_gamma_overtime,dhip_gamma_overtime,'+r')
            hold on
        end
        plot(pl_gamma_overtime_1,dhip_gamma_overtime_2,'k+');
        title('PL-dHPC gamma power correlation')
        xlabel('PL gamma power (mV^2/Hz)')
        ylabel('dHPC gamma power (mV^2/Hz)')
        box off

        subplot(2,3,6) %IL-PL
        if outopt
            plot(il_gamma_overtime,pl_gamma_overtime,'+r')
            hold on
        end
        plot(il_gamma_overtime_3,pl_gamma_overtime_3,'k+');
        title('IL-PL gamma power correlation')
        xlabel('IL gamma power (mV^2/Hz)')
        ylabel('PL gamma power (mV^2/Hz)')
        box off
        
        % Create textboxes displaying Rsquared values (outliers removed)
        %IL-VHIP gamma
        annotation(gcf,'textbox','String',['r^2= ', num2str(Rsquare.Gamma.(stage).il_vhip(idx))],'FitHeightToText','off',...
            'LineStyle','none','Position',[0.2711 0.5888 0.0650 0.0345]);

        % IL-DHIP  gamma
        annotation(gcf,'textbox','String',['r^2= ', num2str(Rsquare.Gamma.(stage).il_dhip(idx))],'FitHeightToText','off',...
            'LineStyle','none','Position',[0.5538 0.5900 0.0629 0.0345]);
        % DHIP-VHIP gamma
        annotation(gcf,'textbox','String',['r^2= ', num2str(Rsquare.Gamma.(stage).dhip_vhip(idx))],'FitHeightToText','off',...
            'LineStyle','none',...
            'Position', [0.8441 0.5875 0.0684 0.0319]);
        % PL-VHIP gamma
        annotation(gcf,'textbox','String',['r^2= ', num2str(Rsquare.Gamma.(stage).pl_vhip(idx))],'FitHeightToText','off',...
            'LineStyle','none',...
            'Position',[0.2649 0.1149 0.0712 0.0370]);
        % PL-DHIP gamma
        annotation(gcf,'textbox','String',['r^2= ', num2str(Rsquare.Gamma.(stage).pl_dhip(idx))],'FitHeightToText','off',...
            'LineStyle','none',...
            'Position',[0.5448 0.1175 0.0698 0.0307]);
        % IL-PL gamma
        annotation(gcf,'textbox','String',['r^2= ', num2str(Rsquare.Gamma.(stage).il_pl(idx))],'FitHeightToText','off',...
            'LineStyle','none',...
            'Position',[0.8232 0.1124 0.0747 0.0409]);   
end %function