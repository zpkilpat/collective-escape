function apply_axes(ax, xl, yl)
% APPLY_AXES  Shared axis cosmetics.
st = paper_style();
xlabel(ax, xl, 'Interpreter','latex', 'FontSize', 36);
ylabel(ax, yl, 'Interpreter','latex', 'FontSize', 36);
set(ax, 'FontSize', st.FSZ, 'TickDir','out', 'LineWidth', 2, ...
        'TickLabelInterpreter','latex', 'Box','off');
end

