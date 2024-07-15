# (1)コンパイラ
CC  = gcc

# (2)コンパイルオプション
CFLAGS		=  -Iinc

# (3)実行ファイル名
TARGET  = calc_Mutual_Information_poisson_noise

# (4)コンパイル対象のソースコード
SRCS    = calc_Mutual_Information_poisson_noise.c inc/binning.c inc/poisson_noise.c

# (5)オブジェクトファイル名
OBJS    = $(SRCS:.c=.o)

# (6)インクルードファイルのあるディレクトリパス
INCDIR  = -I../inc

# (7)ライブラリファイルのあるディレクトリパス
LIBDIR  = 

# (8)追加するライブラリファイル
LIBS    = 

# (9)ターゲットファイル生成
$(TARGET): $(OBJS)
	$(CC) -o $@ $^ $(LIBDIR) $(LIBS)
	@rm -f $(OBJS)

# (10)オブジェクトファイル生成
%.o: %.c
	$(CC) $(CFLAGS) $(INCDIR) -c $< -o $@

# (11)"make all"で make cleanとmakeを同時に実施。
all: $(TARGET) clean_obj

# (12).oファイル、実行ファイル、.dファイルを削除
clean:
	-rm -f $(OBJS) $(TARGET) *.d

# (13) .oファイルのみ削除
clean_obj:
	-rm -f $(OBJS)
