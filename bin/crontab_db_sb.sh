
# test
#*/2 * * * * echo "2" >> /tmp/test2.log

# every six months updates the databases
0 0 1 */6 * "d:/projects/databases/bin/create_db_sb.sh" &>> "d:/projects/databases/logs/create_db_sb.log"
