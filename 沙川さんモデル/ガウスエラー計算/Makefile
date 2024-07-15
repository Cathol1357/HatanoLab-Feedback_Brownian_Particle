# (1)コンパイラ
CC  = gcc

# (2)コンパイルオプション
CFLAGS		=  -Wall -g -Iinc

# (3)実行ファイル名
TARGET  = Feedback_gaussError

# (4)コンパイル対象のソースコード
SRCS    = Feedback_gaussError.c inc/binning.c 

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
	
# (10)オブジェクトファイル生成
%.o: %.c
	$(CC) $(CFLAGS) $(INCDIR) -c $< -o $@

# (11)"make all"で make cleanとmakeを同時に実施。
all: clean $(TARGET)

# (12).oファイル、実行ファイル、.dファイルを削除
clean:
	-rm -f $(OBJS) $(TARGET) *.d
