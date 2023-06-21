function[] = setbgcolor(hand,farb)
set(hand,'BackgroundColor',[0.831 0.816 0.784]);
aa=get(hand,'Units');
set(hand,'Units','pixels');
a=int16(get(hand,'position'));
d(:,:,1)=ones(a(4)-4,a(3)-4)*farb(1);
d(:,:,2)=ones(a(4)-4,a(3)-4)*farb(2);
d(:,:,3)=ones(a(4)-4,a(3)-4)*farb(3);
set(hand,'CData',d);
set(hand,'Units',aa);