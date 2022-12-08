function back = save2txt( file_Name, matrix ) %把矩阵matrix保存成txt文件。

fop = fopen( file_Name, 'wt' );
[M,N] = size(matrix);
for n = 1:N
    for m = 1:M
fprintf( fop, ' %s', mat2str( matrix(m,n) ) );
    end
fprintf(fop, '\n' );
end
back = fclose( fop ) ;