function back = save2txt_1( file_Name, matrix ) %�Ѿ���matrix�����txt�ļ���

fop = fopen( file_Name, 'wt' );
[M,N,A] = size(matrix);

for a=1:A
        for n = 1:N
            for m = 1:M
                fprintf( fop, ' %s', mat2str( matrix(m,n,a) ) );
            end
            fprintf(fop, '\n' );
        end
end