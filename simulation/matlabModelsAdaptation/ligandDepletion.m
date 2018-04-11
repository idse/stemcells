figure; hold on;
Ldose = 0:2:20;
for jj=1:length(Ldose)
    ii = Ldose(jj);
    sol = ode23(@ligandDepletionder,[0 20],[0 1 ii]);
    plot(sol.x,sol.y(1,:),'.-','LineWidth',3);
    leg{jj}=int2str(ii);
end
legend(leg,'Location','Best');
xlabel('time','FontSize',24);
ylabel('Intensity','FontSize',24);
export_fig('ligandDepletion.eps','-transparent');

function dx = ligandDepletionder(t,x)
g = 1; d = 1; b = 1; di = 0.1;
dx = zeros(3,1);

dx(1) = x(3)*x(2) -(g+d)*x(1);
dx(2) = b - (di+x(3))*x(2)+g*x(1);
dx(3) = -x(3)*x(2)+g*x(1);
end