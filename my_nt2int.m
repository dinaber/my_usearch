function SEQ = my_nt2int(NT)
%             A B C D E F G H I J K L M N O P Q R S T U V W X Y Z
nt2i = uint8([1 0 2 0 0 0 3 1 0 0 3 0 1 3 0 0 0 1 2 4 4 0 1 0 2 0 ]) ;
    
SEQ = nt2i(NT-65+1) ;

end
