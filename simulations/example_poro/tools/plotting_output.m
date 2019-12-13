% version='$Rev: 20 $ ($Date: 2017-06-07 18:22:15 +0200 (Mi, 07 Jun 2017) $, $Author: Janis Heuel, Marc S. Boxberg $)'
%% This script plots the output of calculation of wavespeed and inverse quality factor.

%##########################################################
% YOU CAN CHANGE FOR OTHER AXIS LABELS AND UNITS THE FOLLOWING LINES
xlabel_unit = Frequency;          % unit of x-axis
xaxis_label = 'Frequency f [Hz]'; % label of x-axis

%xlabel_unit = omega;
%xaxis_label = 'Angular frequency \omega [rad/s]';

y_axis_label_wavespeed = 'Wavespeed [m/s]';
y_axis_label_Q         = 'Q^{-1}';
%##########################################################

figure(figure_number)
set(figure(figure_number),'name',['Wavespeeds of material ' num2str(u) ' and source ' num2str(src)])
if attenuation == true
    figure(figure_number+1)
    set(figure(figure_number+1),'name',['Attenuation of material ' num2str(u) ' and source ' num2str(src)])
end
if number_of_fluids == 0 || phi_in == 0
    for i = 1:2
        figure(figure_number)
        subplot(2,1,i)
        if i == 1
            semilogx(xlabel_unit,v1)
            legend('P')
        elseif i == 2
            semilogx(xlabel_unit,v2)
            legend('S')
        end
        ylabel(y_axis_label_wavespeed)
        xlabel(xaxis_label)
        if save_plot == true
            print(['w' num2str(u) 's' num2str(src)],formattype)
        end
        if attenuation == true
            figure(figure_number+1)
            subplot(2,1,i)
            if i == 1
                semilogx(xlabel_unit,Q1)
                legend('P')
            elseif i == 2
                semilogx(xlabel_unit,Q2)
                legend('S')
            end
            ylabel(y_axis_label_Q)
            xlabel(xaxis_label)
            if save_plot == true
                print(['a' num2str(u) 's' num2str(src)],formattype)
            end
        end
    end
elseif number_of_fluids == 1
    for i = 1:3
        figure(figure_number)
        subplot(3,1,i)
        if i == 1
            semilogx(xlabel_unit,v1)
            legend('P1')
        elseif i == 2
            semilogx(xlabel_unit,v2)
            legend('S')
        else
            semilogx(xlabel_unit,v3)
            legend('P2')
        end
        ylabel(y_axis_label_wavespeed)
        xlabel(xaxis_label)
        if save_plot == true
            print(['w' num2str(u) 's' num2str(src)],formattype)
        end
        if attenuation == true
            figure(figure_number+1)
            subplot(3,1,i)
            if i == 1
                semilogx(xlabel_unit,Q1)
                legend('P1')
            elseif i == 2
                semilogx(xlabel_unit,Q2)
                legend('S')
            elseif i == 3
                semilogx(xlabel_unit,Q3)
                legend('P2')
            end
            ylabel(y_axis_label_Q)
            xlabel(xaxis_label)
            if save_plot == true
                print(['a' num2str(u) 's' num2str(src)],formattype)
            end
        end
    end
else
    for i = 1:4
        figure(figure_number)
        subplot(2,2,i)
        if i == 1
            semilogx(xlabel_unit,v1)
            legend('P1')
        elseif i == 2
            semilogx(xlabel_unit,v2)
            legend('S')
        elseif i == 3
            semilogx(xlabel_unit,v3)
            legend('P2')
        else
            semilogx(xlabel_unit,v4)
            legend('P3')
        end
        xlabel(xaxis_label)
        ylabel(y_axis_label_wavespeed)
        if save_plot == true
            print(['w' num2str(u) 's' num2str(src)],formattype)
        end
        if attenuation == true
            figure(figure_number+1)
            subplot(2,2,i)
            if i == 1
                semilogx(xlabel_unit,Q1)
                legend('P1')
            elseif i == 2
                semilogx(xlabel_unit,Q2)
                legend('S')
            elseif i == 3
                semilogx(xlabel_unit,Q3)
                legend('P2')
            elseif i == 4
                semilogx(xlabel_unit,Q4)
                legend('P3')
            end
            ylabel(y_axis_label_Q)
            xlabel(xaxis_label)
            if save_plot == true
                print(['a' num2str(u) 's' num2str(src)],formattype)
            end
        end
    end
end
