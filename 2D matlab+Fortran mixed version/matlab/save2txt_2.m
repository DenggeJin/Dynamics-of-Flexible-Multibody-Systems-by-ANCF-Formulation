function back = save2txt_2( file_Name, matrix ) %把矩阵matrix保存成txt文件。

fop = fopen( file_Name, 'wt' );
[M,N,A,B] = size(matrix);
%AB顺序为11，12，13....21，22，23......31，32，33.....

for a=1:A
    for b=1:B
        for n = 1:N
            for m = 1:M
                fprintf( fop, ' %s', mat2str( matrix(m,n,a,b) ) );
            end
            fprintf(fop, '\n' );
        end
    end
end

back = fclose( fop ) ;