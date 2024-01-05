import pandas as pd
import rowindex as rw
import sqlite3
from  sqlite3 import Error


def create_connection(database):
    con = None
    try:
        con= sqlite3.connect(database)
    except Error as e:
        print(e)
    return con

def createTable(con):
    sql = '''CREATE TABLE ModelBlock (
                            ModelName,
                            Volumen,
                            Density,
                            Cu,
                            Litho,
                            OreSort,
                            DbType,
                            ProjectModel_ModelName
                        )'''
    cur = con.cursor()
    cur.execute(sql)
    con.commit()

def t_populateEmpty(con, blocks):
    ModelName = "AS_Model"
    EmptyValues = (ModelName,-1,-1,-1," "," ",-1,ModelName)

    sql = '''INSERT INTO ModelBlock (
                            ModelName,
                            Volumen,
                            Density,
                            Cu,
                            Litho,
                            OreSort,
                            DbType,
                            ProjectModel_ModelName
                        )
                        VALUES (?,?,?,?,?,?,?,?);'''
    cur =con.cursor()
    for i in range(blocks):
        cur.execute(sql,EmptyValues)
    con.commit()
    print(cur.lastrowid)
    return cur.lastrowid

def loadTable(con, bm_csv):
    data = pd.read_csv(bm_csv,sep=',',dtype={"X":float,"Y":float,"Z":float,
        "vol":float, "Cu":float, "density":float, "litho":str, "Ore_sort":str})
    print("Table for update: ", data.shape)

    sql =  '''
    UPDATE ModelBlock 
        SET 
        Volumen = ?,
        Cu = ?,
        Density = ?,
        Litho = ?,
        OreSort = ?,
        DbType = ?
    WHERE
        rowid = ?;
    '''
    try:
        cur = con.cursor()
        for row in data.itertuples(index=False):
            point, values = tuple(row)[:3], tuple(row)[3:]
            index_class = rw.Rowindex(point)
            result = values + (1, index_class.idRow())
            # print(point, result)
            cur.execute(sql,result)
        con.commit()
    except Error as error:
        print("Failed to update Modelblocks", error)
    finally:
        print("The table was updated")
    

def main():
    database = r"/root/kaz/kaz-simulate-grade/extract_weights/data/modelblocks.db"
    bm_csv = r"/root/kaz/kaz-simulate-grade/extract_weights/data/Bm.csv"
    con = create_connection(database)

    with con:
        numberblocks = 4750000
        # Create Table
        # createTable(con)
        # Populate db_empty crea una tabla desde cero
        # t_populateEmpty(con, numberblocks)
        #Actualiza el modelo de bloques de corto plazo
        loadTable(con, bm_csv)

if __name__== '__main__':
    main()