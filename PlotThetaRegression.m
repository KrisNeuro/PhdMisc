function [ThetaPowerCorr] = PlotThetaRegression(outopt,Rsquare,stage,idx, ...
            il_theta_overtime,il_theta_overtime_1,il_theta_overtime_2,il_theta_overtime_3, ...
            dhip_theta_overtime,dhip_theta_overtime_1,dhip_theta_overtime_2,dhip_theta_overtime_3, ...
            vhip_theta_overtime,vhip_theta_overtime_1,vhip_theta_overtime_2,vhip_theta_overtime_3, ...
            pl_theta_overtime,pl_theta_overtime_1,pl_theta_overtime_2,pl_theta_overtime_3)

% Plot theta power regression
ThetaPowerCorr = figure('Color',[1 1 1],'units','normalized','outerposition',[0 0 1 1]);
        subplot(2,3,1) %IL-VHIP
        if outopt
            plot(il_theta_overtime,vhip_theta_overtime,'+r')
            hold on
        end
        plot(il_theta_overtime_2,vhip_theta_overtime_1,'+k')
        title('IL-vHPC theta power correlation')
        xlabel('IL theta power (mV^2/Hz)')
        ylabel('vHPC theta power (mV^2/Hz)')
        box off

        subplot(2,3,2) %IL-DHIP
        if outopt
            plot(il_theta_overtime,dhip_theta_overtime,'+r')
            hold on
        end
        plot(il_theta_overtime_1,dhip_theta_overtime_1,'+k')
        title('IL-dHPC theta power correlation ')
        xlabel('IL theta power (mV^2/Hz)')
        ylabel('dHPC theta power (mV^2/Hz)')
        box off

        subplot(2,3,3) %DHIP-VHIP
        if outopt
            plot(dhip_theta_overtime,vhip_theta_overtime,'+r')
            hold on
        end
        plot(dhip_theta_overtime_3,vhip_theta_overtime_3,'+k')
        title('dHPC-vHPC theta power correlation')
        xlabel('dHPC theta power (mV^2/Hz)')
        ylabel('vHPC theta power (mV^2/Hz)')
        box off

        subplot(2,3,4) %PL-VHIP
        if outopt
            plot(pl_theta_overtime,vhip_theta_overtime,'+r')
            hold on
        end
        plot(pl_theta_overtime_2,vhip_theta_overtime_2,'k+');
        title('PL-vHPC theta power correlation')
        xlabel('PL theta power (mV^2/Hz)')
        ylabel('vHPC theta power (mV^2/Hz)')
        box off

        subplot(2,3,5) %PL-DHIP
        if outopt
            plot(pl_theta_overtime,dhip_theta_overtime,'+r')
            hold on
        end
        plot(pl_theta_overtime_1,dhip_theta_overtime_2,'k+');
        title('PL-dHPC theta power correlation')
        xlabel('PL theta power (mV^2/Hz)')
        ylabel('dHPC theta power (mV^2/Hz)')
        box off

        subplot(2,3,6) %IL-PL
        if outopt
            plot(il_theta_overtime,pl_theta_overtime,'+r')
            hold on
        end
        plot(il_theta_overtime_3,pl_theta_overtime_3,'k+');
        title('IL-PL theta power correlation')
        xlabel('IL theta power (mV^2/Hz)')
        ylabel('PL theta power (mV^2/Hz)')
        box off
        
        % Create textboxes displaying Rsquared values (outliers removed)
if strcmp(subjID(1),'E') %female
            % Create textboxes displaying Rsquared values (outliers removed)
            %IL-VHIP theta
            annotation(gcf,'textbox','String',['r^2= ', num2str(Rsquare.Theta.(stage).il_vhip(idx))],'FitHeightToText','off',...
                'LineStyle','none','Position',[0.2711 0.5888 0.0650 0.0345]);

            % IL-DHIP  theta
            annotation(gcf,'textbox','String',['r^2= ', num2str(Rsquare.Theta.(stage).il_dhip(idx))],'FitHeightToText','off',...
                'LineStyle','none','Position',[0.5538 0.5900 0.0629 0.0345]);
            % DHIP-VHIP theta
            annotation(gcf,'textbox','String',['r^2= ', num2str(Rsquare.Theta.(stage).dhip_vhip(idx))],'FitHeightToText','off',...
                'LineStyle','none',...
                'Position', [0.8441 0.5875 0.0684 0.0319]);
            % PL-VHIP theta
            annotation(gcf,'textbox','String',['r^2= ', num2str(Rsquare.Theta.(stage).pl_vhip(idx))],'FitHeightToText','off',...
                'LineStyle','none',...
                'Position',[0.2649 0.1149 0.0712 0.0370]);
            % PL-DHIP theta
            annotation(gcf,'textbox','String',['r^2= ', num2str(Rsquare.Theta.(stage).pl_dhip(idx))],'FitHeightToText','off',...
                'LineStyle','none',...
                'Position',[0.5448 0.1175 0.0698 0.0307]);
            % IL-PL theta
            annotation(gcf,'textbox','String',['r^2= ', num2str(Rsquare.Theta.(stage).il_pl(idx))],'FitHeightToText','off',...
                'LineStyle','none',...
                'Position',[0.8232 0.1124 0.0747 0.0409]);   
        else %male
             % Create textboxes displaying Rsquared values (outliers removed)
            %IL-VHIP theta
            annotation(gcf,'textbox','String',['r^2= ', num2str(Rsquare.Theta.il_vhip(idx))],'FitHeightToText','off',...
                'LineStyle','none','Position',[0.2711 0.5888 0.0650 0.0345]);
            % IL-DHIP  theta
            annotation(gcf,'textbox','String',['r^2= ', num2str(Rsquare.Theta.il_dhip(idx))],'FitHeightToText','off',...
                'LineStyle','none','Position',[0.5538 0.5900 0.0629 0.0345]);
            % DHIP-VHIP theta
            annotation(gcf,'textbox','String',['r^2= ', num2str(Rsquare.Theta.dhip_vhip(idx))],'FitHeightToText','off',...
                'LineStyle','none',...
                'Position', [0.8441 0.5875 0.0684 0.0319]);
            % PL-VHIP theta
            annotation(gcf,'textbox','String',['r^2= ', num2str(Rsquare.Theta.pl_vhip(idx))],'FitHeightToText','off',...
                'LineStyle','none',...
                'Position',[0.2649 0.1149 0.0712 0.0370]);
            % PL-DHIP theta
            annotation(gcf,'textbox','String',['r^2= ', num2str(Rsquare.Theta.pl_dhip(idx))],'FitHeightToText','off',...
                'LineStyle','none',...
                'Position',[0.5448 0.1175 0.0698 0.0307]);
            % IL-PL theta
            annotation(gcf,'textbox','String',['r^2= ', num2str(Rsquare.Theta.il_pl(idx))],'FitHeightToText','off',...
                'LineStyle','none',...
                'Position',[0.8232 0.1124 0.0747 0.0409]);   


end %function

