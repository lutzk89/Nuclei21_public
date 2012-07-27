function out = filter2s( h, im )
[s1,s2]=size(im);
[sh1,sh2]=size(h);
ial=1; iar=sh1; ibl=sh1+1; ibr=s1+sh1; icl=s1+sh1+1; icr=s1+2*sh1;
jal=1; jar=sh2; jbl=sh2+1; jbr=s2+sh2; jcl=s2+sh2+1; jcr=s2+2*sh2;

out=zeros(icr,jcr);
out(ibl:ibr,jbl:jbr)=im;

out(ial:iar, jbl:jbr)=im(1*ones(sh1,1),:);
out(icl:icr, jbl:jbr)=im(s1*ones(sh1,1),:);

out(ibl:ibr, jal:jar)=im(:,1*ones(1,sh2));
out(ibl:ibr, jcl:jcr)=im(:,s2*ones(1,sh2));

out(ial:iar,jal:jar)=(out(ial+sh1:iar+sh1,jal:jar)+out(ial:iar,jal+sh2:jar+sh2))*.5;
out(icl:icr,jal:jar)=(out(icl-sh1:icr-sh1,jal:jar)+out(icl:icr,jal+sh2:jar+sh2))*.5;
out(ial:iar,jcl:jcr)=(out(ial+sh1:iar+sh1,jcl:jcr)+out(ial:iar,jcl-sh2:jcr-sh2))*.5;
out(icl:icr,jcl:jcr)=(out(icl-sh1:icr-sh1,jcl:jcr)+out(icl:icr,jcl-sh2:jcr-sh2))*.5;

%figure(15);imagesc(out);colorbar;
out=filter2(h,out);
%figure(16);imagesc(out);colorbar;
out=out(ibl:ibr,jbl:jbr);
end
