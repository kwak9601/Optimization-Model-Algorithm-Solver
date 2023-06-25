% Lagrange Algorithm (Equality & In-Equality Constraints)
syms x1 x2 lamda1 mu
f=x1^2+x2^2
disp('First Constraint')
hx1=x1+2*x1+x2-1                 % Equality Constraint 1 (Quadratic + Linear is fine ignore error at the end)
disp('Second Constraint')
gx1=x1-x2                             % Inequality Constraint 2 (Linear + Quadratic) Doesn't Matter
fprintf('Lagrange Equations is:')
l=f+lamda1*(hx1)+mu*(gx1)                % Lagrange Function

% Kuhn-Tucker Conditions
disp('Conditions of Kuhn-tucker Points')
fprintf('1) ')
diff_x1=diff(l,x1)==0;
disp(diff_x1)
fprintf('2) ')
diff_x2=diff(l,x2)==0;
disp(diff_x1)
fprintf('3) ')
disp(hx1==0)
fprintf('4) ')
disp(gx1<=0)
fprintf('5) ')
disp(mu*gx1==0)
fprintf('6) ')
disp(mu>=0)

Y=hessian(l,[x1,x2]);

fprintf('<strong>Case 1:</strong> When %s=0 \n',gx1)
for i=1:2
    clear x1 x2 a lamda1 mu
    syms x1 x2 lamda1 mu
    eq1=gx1==0;
    a=solve(eq1,x2);
    eq2=subs(hx1==0,x2,a);
    fprintf('\n Below r is x1.')
    r=double(solve(eq2,x1,"Real",true))
    disp('At:')
    x1=r(i,1)
    fprintf('At x1=%.4f value of x2 is:',x1)
    x2=double(subs(a,x1))
    eq3=subs(diff_x1);
    eq4=subs(diff_x2);
    [e,f]=solve([eq3,eq4],[lamda1,mu]);
    lamda1=double(e)
    mu=double(f)
    if mu<0
        fprintf('As point mu=%.4f <0 \n ',mu)
        fprintf('The Points %.4f %.4f are not Kuhn-Tuckr Points \n',x1,x2)
    else
        fprintf('The Points %.4f %.4f are Kuhn-Tucker, but lets check the SOSC\n',x1,x2)
        clear x1 x2
        syms x1 x2
        fprintf('<strong>Second Order Sufficient Condition</strong>!\n See Below:\n')
        fprintf(['Second Order Sufficient Conditions Says:' ...
        '\n 2) For every non-zero vector y ∈ T~(x*, μ*), we have yTHx(x*, λ*, μ*) y > 0,'])
        clear y1 y2 y t  P
        syms y1 y2
        t=[y1;y2];
        x1=r(i,1)
        x2=double(subs(a,x1))
        p=subs(gradient(hx1));
        Dh_x=double([p'])
        Dh_xy=Dh_x*t   
        y2=solve(Dh_xy,y2)
    if isempty(y2)
        y2=0
    end
        t=[y1;y2]
        H=subs(Y)
        yHy=t'*H*t;
        yHy=simplify(yHy);
        yHy=subs(yHy)
        y1=rand()
        y2=rand()
    if subs(yHy)>0
        disp((subs(yHy))>0)
        disp('The answer is always positive ')
        fprintf('Hence, The Kuhn-Tucker points %.4f %.4f are strict local minimzers',x1,x2)
    else
        disp('It does not satisfy the SOSC')
    end
    end
end

fprintf('<strong>Case 2:</strong> When mu=0 \n')
for i=1:2
    clear x1 x2
    syms x1 x2
    mu=0;
    eq1=diff_x1==diff_x2;
    a=solve([eq1,hx1],[x1,x2]);
    disp('At:')
    x1=subs(a.x1);
    x2=subs(a.x2);
    x1=x1(i,1)
    x2=x2(i,1)
    constraint=subs(gx1)
if constraint<=0
    fprintf('As the constraint g(x1,x2)= %.4f <=0.\n',constraint)
    fprintf('Hence the points x1=%.4f and x2=%.4f are Kuhn-Tuker Points\n',x1,x2)
    logic=1;
    if logic==1;
        fprintf('The Points %.4f %.4f are Kuhn-Tucker, but lets check the SOSC\n',x1,x2)
        clear x1 x2
        syms x1 x2
        fprintf('<strong>Second Order Sufficient Condition</strong>!\n See Below:\n')
        fprintf(['Second Order Sufficient Conditions Says:' ...
        '\n 1) For every non-zero vector y ∈ T~(x*, μ*), we have yTHx(x*, λ*, μ*) y > 0,'])
        clear y1 y2 y t  P
        syms y1 y2
        t=[y1;y2];
        x1=subs(a.x1);
        x2=subs(a.x2);
        x1=x1(i,1)
        x2=x2(i,1)
        p=subs(gradient(hx1));
        Dh_x=double([p']);
        Dh_xy=Dh_x*t    
        y2=solve(Dh_xy,y2)
    if isempty(y2)
        y2=0
    end
        t=[y1;y2]
        H=subs(Y)
        yHy=t'*H*t;
        yHy=simplify(yHy);
        yHy=subs(yHy)
        y1=rand()
        y2=rand()
    if subs(yHy)>0
        disp((subs(yHy))>0)
        disp('The answer is always positive ')
        fprintf('Hence, The Kuhn-Tucker points %.4f %.4f are strict local minimzers',x1,x2)
    else
        disp('It does not satisfy the SOSC')
    end
    end
else
    fprintf('As the constraint g(x1,x2)= %.4f is not <=0.\n',constraint)
    fprintf('Hence the points x1=%.4f and x2=%.4f are not Kuhn-Tuker Points \n',x1,x2)
end
end
