function back = save2txt( file_Name, matrix ) %�Ѿ���matrix�����txt�ļ���

fop = fopen( file_Name, 'wt' );
[M,N] = size(matrix);
for n = 1:N
    for m = 1:M
fprintf( fop, ' %s', mat2str( matrix(m,n) ) );
    end
fprintf(fop, '\n' );
end
back = fclose( fop ) ;