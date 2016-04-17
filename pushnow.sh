#!/bin/bash
#if hash zenity 2>/dev/null; then
if [ ! `echo $DISPLAY` == "" ]; then
    CommitComment=$(zenity --entry --text "Care to comment?" --entry-text "generic new commit");
else
  CommitComment="remote commit thru terminal; too lazy to type";
fi

git add --all
git commit -m "$CommitComment"
git push